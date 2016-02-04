/**
 * \brief Constructor.
 */
template <int dim>
TransientExecutioner<dim>::TransientExecutioner(
  const TransportParameters<dim> & parameters_,
  Triangulation<dim> & triangulation_,
  const Tensor<1, dim> & transport_direction_,
  const FunctionParser<dim> & cross_section_function_,
  FunctionParser<dim> & source_function_,
  Function<dim> & incoming_function_,
  FunctionParser<dim> & initial_conditions_function_,
  const double & domain_volume_,
  PostProcessor<dim> & postprocessor_,
  const bool & source_is_time_dependent_)
  : Executioner<dim>(parameters_,
                     triangulation_,
                     transport_direction_,
                     cross_section_function_,
                     source_function_,
                     incoming_function_,
                     domain_volume_,
                     postprocessor_),
    temporal_discretization(this->parameters.temporal_discretization),
    initial_conditions_function(&initial_conditions_function_),
    source_is_time_dependent(source_is_time_dependent_),
    theta(this->parameters.theta)
{
  // initialize matrices
  consistent_mass_matrix.reinit(this->constrained_sparsity_pattern);
  lumped_mass_matrix.reinit(this->constrained_sparsity_pattern);
  high_order_diffusion_matrix_new.reinit(this->constrained_sparsity_pattern);
  high_order_ss_matrix_new.reinit(this->constrained_sparsity_pattern);

  // compute nominal time step size
  dt_nominal = this->parameters.time_step_size;

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
    this->dof_handler, 1, *(this->incoming_function), boundary_values);
  for (std::map<unsigned int, double>::const_iterator it =
         boundary_values.begin();
       it != boundary_values.end();
       ++it)
    this->new_solution(it->first) = (it->second);
  // impose other constraints such as hanging nodes on initial conditions
  this->constraints.distribute(this->new_solution);

  // output initial conditions if user requested
  if (parameters_.output_initial_solution)
  {
    std::stringstream IC_filename_ss;
    IC_filename_ss << "solution_" << parameters_.problem_id << "_initial";
    postprocessor_.output_solution(
      this->new_solution, this->dof_handler, IC_filename_ss.str());
  }

  // initialize transient solution vectors
  old_solution.reinit(this->n_dofs);
  older_solution.reinit(this->n_dofs);
  oldest_solution.reinit(this->n_dofs);
  old_stage_solution.reinit(this->n_dofs);
  ss_rhs_new.reinit(this->n_dofs);
  tmp_vector.reinit(this->n_dofs);
}

/**
 * \brief Runs transient executioner.
 */
template <int dim>
void TransientExecutioner<dim>::run()
{
  // assemble inviscid steady-state matrix
  this->assembleInviscidSteadyStateMatrix();

  // assemble steady-state rhs
  this->assembleSteadyStateRHS(this->ss_rhs, 0.0);
  if (!source_is_time_dependent)
    ss_rhs_new = this->ss_rhs;

  // assemble mass matrices
  assembleMassMatrices();

  // compute low-order viscosity
  LowOrderViscosity<dim> low_order_viscosity(this->n_cells,
                                             this->dofs_per_cell,
                                             this->dof_handler,
                                             this->constraints,
                                             this->inviscid_ss_matrix,
                                             this->low_order_diffusion_matrix,
                                             this->low_order_ss_matrix);

  // enforce CFL condition on nominal time step size
  dt_nominal = enforceCFLCondition(dt_nominal);

  // create entropy viscosity
  EntropyViscosity<dim> EV(this->fe,
                           this->n_cells,
                           this->dof_handler,
                           this->constraints,
                           this->cell_quadrature,
                           this->face_quadrature,
                           this->transport_direction,
                           *this->cross_section_function,
                           *this->source_function,
                           this->parameters.entropy_string,
                           this->parameters.entropy_derivative_string,
                           this->parameters.entropy_residual_coefficient,
                           this->parameters.jump_coefficient,
                           this->domain_volume,
                           this->parameters.entropy_temporal_discretization,
                           low_order_viscosity,
                           this->inviscid_ss_matrix,
                           this->high_order_diffusion_matrix,
                           this->high_order_ss_matrix);
  // create FCT object
  FCT<dim> fct(this->dof_handler,
               *this->triangulation,
               lumped_mass_matrix,
               consistent_mass_matrix,
               this->linear_solver,
               this->constrained_sparsity_pattern,
               this->dirichlet_nodes,
               this->n_dofs,
               this->dofs_per_cell,
               this->parameters.do_not_limit,
               this->parameters.theta);

  // create objects required by temporal integrator
  std::shared_ptr<SSPRKTimeIntegrator<dim>> ssprk;
  if (temporal_discretization == TemporalDiscretization::ssprk)
  {
    ssprk = std::make_shared<SSPRKTimeIntegrator<dim>>(
      this->parameters.ssprk_method,
      this->n_dofs,
      this->linear_solver,
      this->constrained_sparsity_pattern);
  }
  else if (temporal_discretization == TemporalDiscretization::theta)
  {
    //
  }

  const double t_end = this->parameters.end_time;
  double t_new = 0.0;
  double t_old = 0.0;
  double dt = dt_nominal;
  double dt_old = dt_nominal;
  double dt_older = dt_nominal;

  old_solution = this->new_solution;
  older_solution = this->new_solution;
  oldest_solution = this->new_solution;

  // compute initial entropy viscosity
  EV.recompute_high_order_ss_matrix(old_solution,
                                    older_solution,
                                    oldest_solution,
                                    dt_old,
                                    dt_older,
                                    t_old,
                                    this->high_order_diffusion_matrix,
                                    this->high_order_ss_matrix);

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

    std::cout << "   time step " << n << ": t = " << t_old << "->" << t_new
              << std::endl;

    if (temporal_discretization == TemporalDiscretization::ssprk) // SSPRK
    {
      // initialize SSPRK time step
      ssprk->initialize_time_step(old_solution, dt);

      switch (this->parameters.viscosity_option)
      {
        case 0: // Galerkin scheme
        {
          compute_galerkin_solution_ssprk(*ssprk);
          break;
        }
        case 1: // Low-order scheme
        {
          compute_low_order_solution_ssprk(*ssprk);
          break;
        }
        case 2: // Entropy viscosity scheme
        {
          compute_entropy_viscosity_solution_ssprk(
            *ssprk, EV, dt_old, dt_older, t_old);
          break;
        }
        case 3: // Entropy viscosity FCT scheme
        {
          compute_entropy_viscosity_fct_solution_ssprk(
            *ssprk, fct, EV, dt, dt_old);
          break;
        }
        case 4: // Galerkin FCT scheme
        {
          compute_galerkin_fct_solution_ssprk(*ssprk, fct, dt);
          break;
        }
        default:
        {
          Assert(false, ExcNotImplemented());
          break;
        }
      }
    }
    else if (temporal_discretization == TemporalDiscretization::theta) // Theta
    {
      switch (this->parameters.viscosity_option)
      {
        case 0: // Galerkin scheme
        {
          compute_galerkin_solution_theta(dt, t_new);
          break;
        }
        case 1: // Low-order scheme
        {
          compute_low_order_solution_theta(dt, t_new);
          break;
        }
        case 2: // Entropy viscosity scheme
        {
          compute_entropy_viscosity_solution_theta(EV, dt, dt_old, t_new);
          break;
        }
        case 3: // Entropy viscosity FCT scheme
        {
          // compute entropy viscosity solution
          compute_entropy_viscosity_solution_theta(EV, dt, dt_old, t_new);

          // compute FCT solution
          compute_fct_solution_theta(fct, dt);

          break;
        }
        case 4: // Galerkin FCT scheme
        {
          // compute Galerkin solution
          compute_galerkin_solution_theta(dt, t_new);

          // compute FCT solution
          compute_fct_solution_theta(fct, dt);

          break;
        }
        default:
        {
          Assert(false, ExcNotImplemented());
          break;
        }
      }

      // save old quantities
      this->ss_rhs = ss_rhs_new;
      this->high_order_diffusion_matrix.copy_from(
        high_order_diffusion_matrix_new);
      this->high_order_ss_matrix.copy_from(high_order_ss_matrix_new);
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }

    // update old solution, time, and time step size
    oldest_solution = older_solution;
    older_solution = old_solution;
    old_solution = this->new_solution;
    dt_older = dt_old;
    dt_old = dt;
    t_old = t_new;

    // increment time step index
    n++;
  }

  // evaluate errors for convergence study
  this->postprocessor->evaluate_error(
    this->new_solution, this->dof_handler, *this->triangulation);

  // output grid and solution and print convergence results
  this->postprocessor->output_results(
    this->new_solution, this->dof_handler, *this->triangulation);

  // print final solution if specified
  if (this->parameters.print_solution)
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
void TransientExecutioner<dim>::assembleMassMatrices()
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
 * @param[in] dt_proposed proposed time step size for current time step
 * @return new time step size
 */
template <int dim>
double TransientExecutioner<dim>::enforceCFLCondition(
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

  // if computed CFL number is greater than the set limit, then adjust dt
  double dt;
  if (proposed_CFL > this->parameters.CFL_limit)
    dt = this->parameters.CFL_limit / max_speed_dx;
  else
    dt = dt_proposed;

  // return adjusted time step size
  return dt;
}

/**
 * \brief Takes time step with the Galerkin method using an SSPRK method.
 */
template <int dim>
void TransientExecutioner<dim>::compute_galerkin_solution_ssprk(
  SSPRKTimeIntegrator<dim> & ssprk)
{
  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
  {
    // get stage time
    double t_stage = ssprk.get_stage_time();

    // recompute steady-state rhs if it is time-dependent
    if (source_is_time_dependent)
    {
      this->assembleSteadyStateRHS(this->ss_rhs, t_stage);
    }

    // advance by an SSPRK step
    ssprk.step(
      consistent_mass_matrix, this->inviscid_ss_matrix, this->ss_rhs, true);
  }
  // retrieve the final solution
  ssprk.get_new_solution(this->new_solution);
}

/**
 * \brief Takes time step with the Galerkin method using a theta method.
 *
 * \param[in] dt time step size
 * \param[in] t_new new time
 */
template <int dim>
void TransientExecutioner<dim>::compute_galerkin_solution_theta(
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
  this->linear_solver.solve(
    this->system_matrix, this->new_solution, this->system_rhs, true, t_new);
}

/**
 * Takes time step with the low-order method using an SSPRK method.
 */
template <int dim>
void TransientExecutioner<dim>::compute_low_order_solution_ssprk(
  SSPRKTimeIntegrator<dim> & ssprk)
{
  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
  {
    // get stage time
    double t_stage = ssprk.get_stage_time();

    // recompute steady-state rhs if it is time-dependent
    if (source_is_time_dependent)
    {
      this->assembleSteadyStateRHS(this->ss_rhs, t_stage);
    }

    // advance by an SSPRK step
    ssprk.step(lumped_mass_matrix, this->low_order_ss_matrix, this->ss_rhs, true);
  }
  // retrieve the final solution
  ssprk.get_new_solution(this->new_solution);
}

/**
 * \brief Takes time step with the low-order method using a theta method.
 *
 * \param[in] dt time step size
 * \param[in] t_new new time
 */
template <int dim>
void TransientExecutioner<dim>::compute_low_order_solution_theta(
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
  this->linear_solver.solve(
    this->system_matrix, this->new_solution, this->system_rhs, true, t_new);
}

/**
 * \brief Takes time step with the entropy viscosity method using an SSPRK method.
 */
template <int dim>
void TransientExecutioner<dim>::compute_entropy_viscosity_solution_ssprk(
  SSPRKTimeIntegrator<dim> & ssprk,
  EntropyViscosity<dim> & EV,
  const double & dt_old,
  const double & dt_older,
  const double & t_old)
{
  // compute EV only at beginning of time step
  EV.recompute_high_order_ss_matrix(old_solution,
                                    older_solution,
                                    oldest_solution,
                                    dt_old,
                                    dt_older,
                                    t_old,
                                    this->high_order_diffusion_matrix,
                                    this->high_order_ss_matrix);

  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
  {
    // get stage time
    double t_stage = ssprk.get_stage_time();

    // recompute steady-state rhs if it is time-dependent
    if (source_is_time_dependent)
      this->assembleSteadyStateRHS(this->ss_rhs, t_stage);

    // advance by an SSPRK step
    ssprk.step(
      consistent_mass_matrix, this->high_order_ss_matrix, this->ss_rhs, true);
  }
  // retrieve the final solution
  ssprk.get_new_solution(this->new_solution);
}

/**
 * \brief Takes time step with the entropy viscosity method using a theta method.
 */
template <int dim>
void TransientExecutioner<dim>::compute_entropy_viscosity_solution_theta(
  EntropyViscosity<dim> & EV,
  const double & dt,
  const double & dt_old,
  const double & t_new)
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
    // recompute new high-order steady-state matrix AH^(l+1)
    EV.recompute_high_order_ss_matrix(this->new_solution,
                                      old_solution,
                                      older_solution,
                                      dt,
                                      dt_old,
                                      t_new,
                                      high_order_diffusion_matrix_new,
                                      high_order_ss_matrix_new);

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
 * Takes time step with the entropy viscosity method with FCT.
 */
template <int dim>
void TransientExecutioner<dim>::compute_entropy_viscosity_fct_solution_ssprk(
  SSPRKTimeIntegrator<dim> & ssprk,
  FCT<dim> & fct,
  EntropyViscosity<dim> & EV,
  const double & dt,
  const double & dt_old)
{
  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
  {
    // get stage time
    double t_stage = ssprk.get_stage_time();

    // recompute steady-state rhs if it is time-dependent
    if (source_is_time_dependent)
      this->assembleSteadyStateRHS(this->ss_rhs, t_stage);

    // compute Galerkin solution
    ssprk.step(
      consistent_mass_matrix, this->inviscid_ss_matrix, this->ss_rhs, false);

    // get Galerkin solution
    ssprk.get_intermediate_solution(this->new_solution);

    // get old stage solution
    ssprk.get_stage_solution(i, old_stage_solution);

    // recompute high-order steady-state matrix
    EV.recompute_high_order_ss_matrix(this->new_solution,
                                      old_stage_solution,
                                      old_solution,
                                      dt,
                                      dt_old,
                                      t_stage,
                                      this->high_order_diffusion_matrix,
                                      this->high_order_ss_matrix);

    // advance by an SSPRK step
    ssprk.step(
      consistent_mass_matrix, this->high_order_ss_matrix, this->ss_rhs, false);

    // get old stage solution
    ssprk.get_stage_solution(i, old_stage_solution);

    // get intermediate solution
    ssprk.get_intermediate_solution(this->new_solution);

    // perform FCT
    fct.solve_FCT_system_fe(this->new_solution,
                            old_stage_solution,
                            this->low_order_ss_matrix,
                            this->ss_rhs,
                            dt,
                            this->low_order_diffusion_matrix,
                            this->high_order_diffusion_matrix);

    // set stage solution to be FCT solution for this stage
    ssprk.set_intermediate_solution(this->new_solution);

    // finish computing stage solution
    ssprk.complete_stage_solution();
  }
  // retrieve the final solution
  ssprk.get_new_solution(this->new_solution);
}

/**
 * Takes time step with the Galerkin method with FCT.
 */
template <int dim>
void TransientExecutioner<dim>::compute_galerkin_fct_solution_ssprk(
  SSPRKTimeIntegrator<dim> & ssprk, FCT<dim> & fct, const double & dt)
{
  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
  {
    // get stage time
    double t_stage = ssprk.get_stage_time();

    // recompute steady-state rhs if it is time-dependent
    if (source_is_time_dependent)
      this->assembleSteadyStateRHS(this->ss_rhs, t_stage);

    // compute Galerkin solution
    ssprk.step(
      consistent_mass_matrix, this->inviscid_ss_matrix, this->ss_rhs, false);

    // get Galerkin solution
    ssprk.get_intermediate_solution(this->new_solution);

    // get old stage solution
    ssprk.get_stage_solution(i, old_stage_solution);

    // perform FCT
    fct.solve_FCT_system_fe(this->new_solution,
                            old_stage_solution,
                            this->low_order_ss_matrix,
                            this->ss_rhs,
                            dt,
                            this->low_order_diffusion_matrix,
                            this->high_order_diffusion_matrix);

    // set stage solution to be FCT solution for this stage
    ssprk.set_intermediate_solution(this->new_solution);

    // finish computing stage solution
    ssprk.complete_stage_solution();
  }

  // retrieve the final solution
  ssprk.get_new_solution(this->new_solution);
}

/**
 * \brief Computes the FCT solution using a theta method.
 *
 * \param[in] dt  time step size
 */
template <int dim>
void TransientExecutioner<dim>::compute_fct_solution_theta(FCT<dim> & fct,
                                                           const double & dt)
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
                             dt);

    // compute limited flux bounds
    fct.compute_limited_flux_bounds_theta(this->new_solution,
                                          old_solution,
                                          this->low_order_ss_matrix,
                                          ss_rhs_new,
                                          ss_rhs_old,
                                          dt);

    // compute limited flux sums
    fct.compute_limited_fluxes();

    // create system rhs
    this->system_rhs = 0;
    this->system_rhs.add((1.0 - theta) * dt, ss_rhs_old);
    this->system_rhs.add(theta * dt, ss_rhs_new);
    lumped_mass_matrix.vmult(tmp_vector, old_solution);
    this->system_rhs.add(1.0, tmp_vector);
    this->low_order_ss_matrix.vmult(tmp_vector, old_solution);
    this->system_rhs.add(-(1.0 - theta) * dt, tmp_vector);
    this->system_rhs.add(dt, fct.get_flux_correction_vector());

    // apply Dirichlet BC here
    this->applyDirichletBC(
      this->system_matrix, this->system_rhs, this->new_solution);

    // check convergence and perform update if necessary
    converged =
      this->nonlinear_solver.update(this->system_matrix, this->system_rhs);
  }
}
