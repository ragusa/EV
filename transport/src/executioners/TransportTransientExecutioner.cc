/**
 * \brief Constructor.
 */
template <int dim>
TransportTransientExecutioner<dim>::TransportTransientExecutioner(
  const TransportRunParameters & parameters_,
  TransportProblemParameters<dim> & problem_parameters_,
  Triangulation<dim> & triangulation_,
  PostProcessor<dim> & postprocessor_,
  const double & nominal_dt_)
  : TransportExecutioner<dim>(
      parameters_, problem_parameters_, triangulation_, postprocessor_),
    temporal_discretization(parameters_.temporal_discretization),
    initial_conditions_function(
      &(problem_parameters_.initial_conditions_function)),
    dt_nominal(nominal_dt_),
    source_is_time_dependent(problem_parameters_.source_is_time_dependent),
    theta(parameters_.theta)
{
  // initialize matrices
  consistent_mass_matrix.reinit(this->constrained_sparsity_pattern);
  lumped_mass_matrix.reinit(this->constrained_sparsity_pattern);
  high_order_diffusion_matrix_new.reinit(this->constrained_sparsity_pattern);
  high_order_ss_matrix_new.reinit(this->constrained_sparsity_pattern);

  // compute minimum cell diameter for CFL condition
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  minimum_cell_diameter = cell->diameter();
  for (; cell != endc; ++cell)
    minimum_cell_diameter = std::min(minimum_cell_diameter, cell->diameter());

  // interpolate initial conditions
  initial_conditions_function->set_time(0.0);
  VectorTools::interpolate(
    this->dof_handler, *initial_conditions_function, this->new_solution);

  // impose Dirichlet BC on initial conditions
  std::map<unsigned int, double> boundary_values;
  VectorTools::interpolate_boundary_values(
    this->dof_handler, 0, *(this->incoming_function), boundary_values);
  for (std::map<unsigned int, double>::const_iterator it =
         boundary_values.begin();
       it != boundary_values.end();
       ++it)
    this->new_solution(it->first) = (it->second);

  // impose other constraints such as hanging nodes on initial conditions
  this->constraints.distribute(this->new_solution);

  // output initial conditions
  postprocessor_.output_solution(
    this->new_solution, 0.0, this->dof_handler, "solution_initial");

  // initialize transient solution vectors
  old_solution.reinit(this->n_dofs);
  older_solution.reinit(this->n_dofs);
  old_stage_solution.reinit(this->n_dofs);
  ss_rhs_new.reinit(this->n_dofs);
  tmp_vector.reinit(this->n_dofs);
}

/**
 * \brief Runs transient executioner.
 */
template <int dim>
void TransportTransientExecutioner<dim>::run()
{
  // assemble inviscid steady-state matrix
  this->assembleInviscidSteadyStateMatrix();

  // compute low-order steady-state matrix since it is not time-dependent
  // NOTE: arbitrary values are passed because the viscosity should not depend on
  // these parameters
  this->low_order_viscosity->update(
    this->new_solution, this->new_solution, 1.0, 1);
  this->low_order_diffusion->compute_diffusion_matrix(
    this->new_solution,
    this->low_order_viscosity,
    this->low_order_diffusion_matrix);
  this->low_order_ss_matrix.copy_from(this->inviscid_ss_matrix);
  this->low_order_ss_matrix.add(1.0, this->low_order_diffusion_matrix);

  // compute initial old high-order steady-state matrix
  // NOTE: this is only needed when entropy viscosity theta scheme is used
  this->high_order_viscosity->update(
    this->new_solution, this->new_solution, 1.0, 1);
  this->high_order_diffusion->compute_diffusion_matrix(
    this->new_solution,
    this->high_order_viscosity,
    this->high_order_diffusion_matrix);
  this->high_order_ss_matrix.copy_from(this->inviscid_ss_matrix);
  this->high_order_ss_matrix.add(1.0, this->high_order_diffusion_matrix);

  // assemble steady-state rhs
  this->assembleSteadyStateRHS(this->ss_rhs, 0.0);
  if (!source_is_time_dependent)
    ss_rhs_new = this->ss_rhs;

  // assemble mass matrices
  assembleMassMatrices();

  // create FCT object
  std::shared_ptr<FCT<dim>> fct;
  if (this->parameters.scheme == Scheme::fct)
  {
    fct = std::make_shared<FCT<dim>>(this->dof_handler,
                                     *this->triangulation,
                                     lumped_mass_matrix,
                                     consistent_mass_matrix,
                                     this->linear_solver,
                                     this->constrained_sparsity_pattern,
                                     this->dirichlet_nodes,
                                     this->n_dofs,
                                     this->dofs_per_cell,
                                     false,
                                     this->parameters.include_analytic_bounds,
                                     this->fe,
                                     this->cell_quadrature,
                                     *this->cross_section_function,
                                     *this->source_function,
                                     source_is_time_dependent,
                                     this->parameters.theta);
  }

  // create objects required by temporal integrator
  std::shared_ptr<SSPRKTimeIntegrator<dim>> ssprk;
  if (temporal_discretization == TemporalDiscretizationClassification::ssprk)
  {
    ssprk = std::make_shared<SSPRKTimeIntegrator<dim>>(
      this->parameters.ssprk_discretization,
      this->n_dofs,
      this->linear_solver,
      this->constrained_sparsity_pattern);
  }

  // enforce CFL condition on nominal time step size
  dt_nominal = enforceCFLCondition(dt_nominal);

  // update nominal dt reported on convergence table
  this->postprocessor->log_time_step_size(dt_nominal);

  // TODO: check for default end time
  const double t_end = this->parameters.end_time;
  double t_new = 0.0;
  double t_old = 0.0;
  double dt = dt_nominal;
  double dt_old = dt_nominal;

  old_solution = this->new_solution;
  older_solution = this->new_solution;

  // time loop
  unsigned int n = 1; // time step index
  bool in_transient = true;
  while (in_transient)
  {
    // shorten time step size if new time would overshoot end time
    dt = dt_nominal;
    if (t_old + dt >= t_end)
    {
      dt = t_end - t_old;
      in_transient = false;
    }

    // compute new time
    t_new = t_old + dt;

    this->cout1 << "   time step " << n << ": t = " << t_old << "->" << t_new
                << std::endl;

    if (temporal_discretization == TemporalDiscretizationClassification::ssprk)
    {
      // initialize SSPRK time step
      ssprk->initialize_time_step(old_solution, dt);

      // loop over SSPRK stages to compute new solution
      for (unsigned int i = 0; i < ssprk->n_stages; ++i)
      {
        // get stage time
        const double t_stage = ssprk->get_stage_time();

        // determine old stage solution and dt (needed for entropy viscosity)
        double old_stage_dt;
        if (i == 0)
        {
          old_stage_dt = dt_old;
        }
        else
        {
          ssprk->get_stage_solution(i - 1, older_solution);
          old_stage_dt = dt;
        }

        // get old stage solution
        ssprk->get_stage_solution(i, old_stage_solution);

        // compute steady-state rhs vector
        if (source_is_time_dependent)
          this->assembleSteadyStateRHS(this->ss_rhs, t_stage);

        // advance by an SSPRK step
        switch (this->parameters.scheme)
        {
          case Scheme::low:
            this->low_order_viscosity->update(
              old_stage_solution, older_solution, old_stage_dt, n);
            this->low_order_diffusion->compute_diffusion_matrix(
              old_stage_solution,
              this->low_order_viscosity,
              this->low_order_diffusion_matrix);
            ssprk->step(lumped_mass_matrix,
                        this->inviscid_ss_matrix,
                        this->low_order_diffusion_matrix,
                        this->ss_rhs,
                        true);
            break;
          case Scheme::high:
            this->high_order_viscosity->update(
              old_stage_solution, older_solution, old_stage_dt, n);
            this->high_order_diffusion->compute_diffusion_matrix(
              old_stage_solution,
              this->high_order_viscosity,
              this->high_order_diffusion_matrix);
            ssprk->step(consistent_mass_matrix,
                        this->inviscid_ss_matrix,
                        this->high_order_diffusion_matrix,
                        this->ss_rhs,
                        true);
            break;
          case Scheme::fct:
            perform_fct_ssprk_step(t_stage, dt, old_stage_dt, n, fct, *ssprk);
            break;
          default:
            throw ExcNotImplemented();
            break;
        }
      }

      // retrieve the final solution
      ssprk->get_new_solution(this->new_solution);
    }
    else if (temporal_discretization ==
             TemporalDiscretizationClassification::theta)
    {
      switch (this->parameters.scheme)
      {
        case Scheme::low:
          compute_low_order_solution_theta(dt, t_new);
          break;
        case Scheme::high:
          compute_high_order_solution_theta(dt, t_new, n);
          break;
        case Scheme::fct:
          compute_high_order_solution_theta(dt, t_new, n);
          compute_fct_solution_theta(*fct, dt, t_old);
          break;
        default:
          throw ExcNotImplemented();
          break;
      }

      // save old quantities
      this->ss_rhs = ss_rhs_new;
      this->high_order_diffusion_matrix.copy_from(
        high_order_diffusion_matrix_new);
      this->high_order_ss_matrix.copy_from(high_order_ss_matrix_new);
    }
    else
    {
      throw ExcNotImplemented();
    }

    // increment transient counter for post-processor
    this->postprocessor->increment_transient_counter();

    // update old solution, time, and time step size
    older_solution = old_solution;
    old_solution = this->new_solution;
    dt_old = dt;
    t_old = t_new;

    // increment time step index
    n++;
  }

  // evaluate errors for convergence study
  this->postprocessor->evaluate_error(
    this->new_solution, this->dof_handler, *this->triangulation);

  // output grid and solution and print convergence results if in last cycle
  this->postprocessor->output_results_if_last_cycle(
    this->new_solution, this->dof_handler, *this->triangulation);

  // output FCT bounds if requested
  if (this->parameters.output_fct_bounds)
    fct->output_bounds(*(this->postprocessor));

  // print final solution if specified
  if (this->parameters.print_final_solution)
  {
    // set precision and format
    std::cout.precision(10);
    std::cout.setf(std::ios::scientific);

    // print each value of solution
    for (unsigned int j = 0; j < this->n_dofs; ++j)
      std::cout << this->new_solution[j] << std::endl;
  }
}

/**
 * Assembles the consistent and lumped mass matrices.
 */
template <int dim>
void TransportTransientExecutioner<dim>::assembleMassMatrices()
{
  FEValues<dim> fe_values(
    this->fe, this->cell_quadrature, update_values | update_JxW_values);

  std::vector<types::global_dof_index> local_dof_indices(this->dofs_per_cell);

  FullMatrix<double> local_mass_consistent(this->dofs_per_cell,
                                           this->dofs_per_cell);
  FullMatrix<double> local_mass_lumped(this->dofs_per_cell, this->dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    fe_values.reinit(cell);
    cell->get_dof_indices(local_dof_indices);

    local_mass_consistent = 0.0;
    local_mass_lumped = 0.0;

    // compute local contribution
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
        for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
        {
          local_mass_consistent(i, j) += fe_values.shape_value(i, q) *
            fe_values.shape_value(j, q) * fe_values.JxW(q);
          local_mass_lumped(i, i) += fe_values.shape_value(i, q) *
            fe_values.shape_value(j, q) * fe_values.JxW(q);
        }

    // add to global mass matrices with contraints
    this->constraints.distribute_local_to_global(
      local_mass_consistent, local_dof_indices, consistent_mass_matrix);
    this->constraints.distribute_local_to_global(
      local_mass_lumped, local_dof_indices, lumped_mass_matrix);
  }
}

/**
 * Checks CFL condition and adjusts time step size as necessary.
 *
 * \param[in] dt_proposed proposed time step size for current time step
 *
 * \return new time step size
 */
template <int dim>
double TransportTransientExecutioner<dim>::enforceCFLCondition(
  const double & dt_proposed) const
{
  // CFL is dt*speed/dx
  double max_speed_dx = 0.0; // max(speed/dx); max over all i of A(i,i)/mL(i,i)
  for (unsigned int i = 0; i < this->n_dofs; ++i)
    max_speed_dx = std::max(max_speed_dx,
                            std::abs(this->low_order_ss_matrix(i, i)) /
                              lumped_mass_matrix(i, i));

  // compute proposed CFL number
  double proposed_CFL = dt_proposed * max_speed_dx;

  double dt;
  // if time step violates CFL limit, then respond accordingly
  if (proposed_CFL > this->parameters.cfl)
    // if user specified to adjust dt based on CFL
    if (this->parameters.time_step_size_option ==
        TimeStepSizeOption::cfl)
    {
      // adjust dt to satisfy CFL
      dt = this->parameters.cfl / max_speed_dx;
    }
    else
    {
      // report CFL violation but do not change dt
      std::cout << "CFL violated for cycle: ";
      std::cout << "\x1b[1;31m CFL = " << proposed_CFL << "\x1b[0m" << std::endl;
      dt = dt_proposed;
    }
  else
    dt = dt_proposed;

  // return adjusted time step size
  return dt;
}

/**
 * \brief Takes time step with the Galerkin method using a theta method.
 *
 * \param[in] dt time step size
 * \param[in] t_new new time
 */
template <int dim>
void TransportTransientExecutioner<dim>::compute_galerkin_solution_theta(
  const double & dt, const double & t_new)
{
  // compute new steady-state right-hand-side vector b_new
  if (source_is_time_dependent)
    this->assembleSteadyStateRHS(ss_rhs_new, t_new);

  // compute system matrix
  // matrix = M^C + theta*dt*A
  this->system_matrix.copy_from(consistent_mass_matrix);
  this->system_matrix.add(theta * dt, this->inviscid_ss_matrix);

  // compute system right-hand-side vector
  // rhs = M^C*u_old + dt*((1-theta)*b_old + theta*b_new - (1-theta)*A*u_old)
  consistent_mass_matrix.vmult(this->system_rhs, old_solution);
  this->inviscid_ss_matrix.vmult(tmp_vector, old_solution);
  this->system_rhs.add(-(1.0 - theta) * dt, tmp_vector);
  this->system_rhs.add((1.0 - theta) * dt, this->ss_rhs);
  this->system_rhs.add(theta * dt, ss_rhs_new);

  // solve linear system
  this->linear_solver.solve_with_dirichlet(
    this->system_matrix, this->new_solution, this->system_rhs, true, t_new);
}

/**
 * \brief Takes time step with the low-order method using a theta method.
 *
 * \param[in] dt     new time step size \f$\Delta t^n\equiv t^{n+1}-t^n\f$
 * \param[in] t_new  new time \f$t^{n+1}\f$
 */
template <int dim>
void TransportTransientExecutioner<dim>::compute_low_order_solution_theta(
  const double & dt, const double & t_new)
{
  // compute new steady-state right-hand-side vector b_new
  if (source_is_time_dependent)
    this->assembleSteadyStateRHS(ss_rhs_new, t_new);

  // compute system matrix
  // matrix = M^L + theta*dt*A^L
  this->system_matrix.copy_from(lumped_mass_matrix);
  this->system_matrix.add(theta * dt, this->low_order_ss_matrix);

  // compute system right-hand-side vector
  // rhs = M^L*u_old + dt*((1-theta)*b_old + theta*b_new - (1-theta)*A^L*u_old)
  lumped_mass_matrix.vmult(this->system_rhs, old_solution);
  this->low_order_ss_matrix.vmult(tmp_vector, old_solution);
  this->system_rhs.add(-(1.0 - theta) * dt, tmp_vector);
  this->system_rhs.add((1.0 - theta) * dt, this->ss_rhs);
  this->system_rhs.add(theta * dt, ss_rhs_new);

  // solve linear system
  this->linear_solver.solve_with_dirichlet(
    this->system_matrix, this->new_solution, this->system_rhs, true, t_new);
}

/**
 * \brief Takes time step with a high-order method using a theta method.
 *
 * \param[in] dt     new time step size \f$\Delta t^n\equiv t^{n+1}-t^n\f$
 * \param[in] t_new  new time \f$t^{n+1}\f$
 * \param[in] n      time step index
 */
template <int dim>
void TransportTransientExecutioner<dim>::compute_high_order_solution_theta(
  const double & dt, const double & t_new, const unsigned int & n)
{
  switch (this->parameters.high_order_scheme)
  {
    case HighOrderScheme::galerkin: // galerkin
      compute_galerkin_solution_theta(dt, t_new);
      break;
    case HighOrderScheme::entropy_visc: // entropy viscosity
      compute_entropy_viscosity_solution_theta(dt, t_new, n);
      break;
    default:
      throw ExcNotImplemented();
      break;
  }
}

/**
 * \brief Takes time step with the entropy viscosity method using a theta method.
 *
 * \param[in] dt     new time step size \f$\Delta t^n    \equiv t^{n+1}-t^n\f$
 * \param[in] t_new  new time \f$t^{n+1}\f$
 * \param[in] n      time step index
 */
template <int dim>
void TransportTransientExecutioner<dim>::compute_entropy_viscosity_solution_theta(
  const double & dt, const double & t_new, const unsigned int & n)
{
  // references
  Vector<double> & ss_rhs_old = this->ss_rhs;
  SparseMatrix<double> & high_order_ss_matrix_old = this->high_order_ss_matrix;

  // compute new steady-state right-hand-side vector b_new
  if (source_is_time_dependent)
    this->assembleSteadyStateRHS(ss_rhs_new, t_new);

  // initialize guess for nonlinear solver
  this->new_solution = 0.0;
  this->nonlinear_solver.initialize(this->new_solution);

  // begin iteration
  bool converged = false;
  while (!converged)
  {
    // compute new entropy viscosity and high order diffusion matrix
    this->high_order_viscosity->update(
      this->new_solution, old_solution, dt, n + 1);
    this->high_order_diffusion->compute_diffusion_matrix(
      this->new_solution,
      this->high_order_viscosity,
      high_order_diffusion_matrix_new);

    // compute new high-order steady-state matrix
    high_order_ss_matrix_new.copy_from(this->inviscid_ss_matrix);
    high_order_ss_matrix_new.add(1.0, high_order_diffusion_matrix_new);

    // compute system matrix
    // matrix = M^C + theta*dt*A^H
    this->system_matrix.copy_from(consistent_mass_matrix);
    this->system_matrix.add(theta * dt, high_order_ss_matrix_new);

    // compute system right-hand-side vector
    // rhs = M^C*u_old + dt*((1-theta)*b_old + theta*b_new - (1-theta)*A^H*u_old)
    consistent_mass_matrix.vmult(this->system_rhs, old_solution);
    high_order_ss_matrix_old.vmult(tmp_vector, old_solution);
    this->system_rhs.add(-(1.0 - theta) * dt, tmp_vector);
    this->system_rhs.add((1.0 - theta) * dt, ss_rhs_old);
    this->system_rhs.add(theta * dt, ss_rhs_new);

    // apply Dirichlet BC
    this->applyDirichletBC(
      this->system_matrix, this->system_rhs, this->new_solution);

    // check convergence and perform update if necessary
    converged =
      this->nonlinear_solver.update(this->system_matrix, this->system_rhs);
  }
}

/**
 * \brief Performs an SSPRK step using FCT.
 *
 * \param[in] t_old old time
 * \param[in] dt current time step size
 * \param[in] old_stage_dt time step size of previous SSPRK stage
 * \param[in] n time index
 * \param[in] fct FCT
 * \param[in,out] ssprk SSPRK time integrator
 */
template <int dim>
void TransportTransientExecutioner<dim>::perform_fct_ssprk_step(
  const double & t_old,
  const double & dt,
  const double & old_stage_dt,
  const unsigned int & n,
  const std::shared_ptr<FCT<dim>> & fct,
  SSPRKTimeIntegrator<dim> & ssprk)
{
  // update low-order diffusion
  this->low_order_viscosity->update(
    old_stage_solution, older_solution, old_stage_dt, n);
  this->low_order_diffusion->compute_diffusion_matrix(
    old_stage_solution,
    this->low_order_viscosity,
    this->low_order_diffusion_matrix);

  // update high-order diffusion
  this->high_order_viscosity->update(
    old_stage_solution, older_solution, old_stage_dt, n);
  this->high_order_diffusion->compute_diffusion_matrix(
    old_stage_solution,
    this->high_order_viscosity,
    this->high_order_diffusion_matrix);

  // compute high-order solution
  ssprk.step(consistent_mass_matrix,
             this->inviscid_ss_matrix,
             this->high_order_diffusion_matrix,
             this->ss_rhs,
             false);

  // get high-order solution
  ssprk.get_intermediate_solution(this->new_solution);

  // perform FCT
  fct->solve_FCT_system_fe(this->new_solution,
                           old_stage_solution,
                           this->low_order_ss_matrix,
                           this->ss_rhs,
                           dt,
                           this->low_order_diffusion_matrix,
                           this->high_order_diffusion_matrix,
                           t_old);

  // set stage solution to be FCT solution for this stage
  ssprk.set_intermediate_solution(this->new_solution);

  // finish computing stage solution
  ssprk.complete_stage_solution();
}

/**
 * \brief Computes the FCT solution using a theta method.
 *
 * \param[in] dt  time step size
 */
template <int dim>
void TransportTransientExecutioner<dim>::compute_fct_solution_theta(
  FCT<dim> & fct, const double & dt, const double & t_old)
{
  // references
  const Vector<double> & ss_rhs_old = this->ss_rhs;
  const SparseMatrix<double> & high_order_diffusion_matrix_old =
    this->high_order_diffusion_matrix;

  // compute flux corrections
  fct.compute_flux_corrections_theta(this->new_solution,
                                     old_solution,
                                     dt,
                                     this->low_order_diffusion_matrix,
                                     high_order_diffusion_matrix_old,
                                     high_order_diffusion_matrix_new);

  // compute system matrix
  this->system_matrix.copy_from(lumped_mass_matrix);
  this->system_matrix.add(theta * dt, this->low_order_ss_matrix);

  // initialize guess for nonlinear solver
  this->new_solution = 0.0;
  this->nonlinear_solver.initialize(this->new_solution);

  // initialize cumulative antidiffusion vector
  this->cumulative_antidiffusion = 0.0;

  // begin iteration
  bool converged = false;
  while (!converged)
  {
    // compute max principle min and max values
    fct.compute_bounds_theta(this->new_solution,
                             old_solution,
                             this->low_order_ss_matrix,
                             ss_rhs_new,
                             ss_rhs_old,
                             dt,
                             t_old);

    // compute limited flux bounds
    fct.compute_limited_flux_bounds_theta(this->new_solution,
                                          old_solution,
                                          this->low_order_ss_matrix,
                                          ss_rhs_new,
                                          ss_rhs_old,
                                          this->cumulative_antidiffusion,
                                          dt);

    // compute limited flux sums
    fct.compute_limited_fluxes();

    // if using cumulative antidiffusion algorithm, then update cumulative
    // antidiffusion and remainder antidiffusive fluxes
    if (this->parameters.use_cumulative_antidiffusion_algorithm)
    {
      // add to cumulative antidiffusion vector: p^(l+1) = p^(l) + dp^(l)
      this->cumulative_antidiffusion.add(1.0, fct.get_limited_flux_vector());

      // subtract used antidiffusive flux: DP^(l+1) = DP^(l) - dP^(l)
      fct.subtract_limited_flux_correction_matrix();
    }
    else
    {
      // throw away accumulated antidiffusive flux
      this->cumulative_antidiffusion = fct.get_limited_flux_vector();
    }

    // create system rhs
    this->system_rhs = 0;
    this->system_rhs.add((1.0 - theta) * dt, ss_rhs_old);
    this->system_rhs.add(theta * dt, ss_rhs_new);
    lumped_mass_matrix.vmult(tmp_vector, old_solution);
    this->system_rhs.add(1.0, tmp_vector);
    this->low_order_ss_matrix.vmult(tmp_vector, old_solution);
    this->system_rhs.add(-(1.0 - theta) * dt, tmp_vector);
    this->system_rhs.add(dt, this->cumulative_antidiffusion);

    // create system matrix. Note that although the underlying system matrix
    // does not change in each iteration, the Dirichlet-modified system matrix
    // does change
    this->system_matrix.copy_from(lumped_mass_matrix);
    this->system_matrix.add(theta * dt, this->low_order_ss_matrix);

    // apply Dirichlet BC here
    this->applyDirichletBC(
      this->system_matrix, this->system_rhs, this->new_solution);

    // check convergence and perform update if necessary
    converged =
      this->nonlinear_solver.update(this->system_matrix, this->system_rhs);
  }

  // check FCT solution
  fct.check_fct_bounds(this->new_solution);
}
