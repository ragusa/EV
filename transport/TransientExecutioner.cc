/**
 * Constructor.
 */
template<int dim>
TransientExecutioner<dim>::TransientExecutioner(
  const TransportParameters<dim> & parameters_,
  const Triangulation<dim> & triangulation_,
  const Tensor<1, dim> & transport_direction_,
  const FunctionParser<dim> & cross_section_function_,
  FunctionParser<dim> & source_function_,
  Function<dim> & incoming_function_,
  FunctionParser<dim> & initial_conditions_function_,
  PostProcessor<dim> & postprocessor_,
  const bool & source_is_time_dependent_) :
    Executioner<dim>(parameters_, triangulation_, transport_direction_,
      cross_section_function_, source_function_, incoming_function_, postprocessor_),
    initial_conditions_function(& initial_conditions_function_),
    source_is_time_dependent(source_is_time_dependent_)
{
  // initialize mass matrices
  consistent_mass_matrix.reinit(this->constrained_sparsity_pattern);
  lumped_mass_matrix.reinit(this->constrained_sparsity_pattern);

  // compute nominal time step size
  dt_nominal = this->parameters.time_step_size;

  // compute minimum cell diameter for CFL condition
  typename DoFHandler<dim>::active_cell_iterator cell =
      this->dof_handler.begin_active(), endc = this->dof_handler.end();
  minimum_cell_diameter = cell->diameter();
  for (; cell != endc; ++cell)
    minimum_cell_diameter = std::min(minimum_cell_diameter, cell->diameter());

  // interpolate initial conditions
  initial_conditions_function->set_time(0.0);
  VectorTools::interpolate(this->dof_handler, *initial_conditions_function, this->new_solution);
  this->constraints.distribute(this->new_solution);

  // if last cycle, output initial conditions if user requested
  /*
  if ((cycle == parameters.n_refinement_cycles - 1)
      and (parameters.output_initial_solution))
  {
    std::stringstream IC_filename_ss;
    IC_filename_ss << "solution_" << parameters.problem_id << "_initial";
    postprocessor.output_solution(new_solution, dof_handler,
        IC_filename_ss.str());
  }
  */

  // initialize transient solution vectors
  old_solution.reinit(this->n_dofs);
  older_solution.reinit(this->n_dofs);
  oldest_solution.reinit(this->n_dofs);
  old_stage_solution.reinit(this->n_dofs);
}

/**
 * Destructor.
 */
template<int dim>
TransientExecutioner<dim>::~TransientExecutioner()
{
}

/**
 * Runs transient executioner.
 */
template<int dim>
void TransientExecutioner<dim>::run()
{
  // assemble inviscid steady-state matrix and steady-state rhs
  this->assembleInviscidSteadyStateMatrix();
  this->assembleSteadyStateRHS(0.0);

  // assemble mass matrices
  assembleMassMatrices();

  // compute low-order viscosity
  LowOrderViscosity<dim> low_order_viscosity(this->n_cells,
    this->dofs_per_cell, this->dof_handler, this->constraints,
    this->inviscid_ss_matrix, this->low_order_diffusion_matrix,
    this->low_order_ss_matrix);

  // enforce CFL condition on nominal time step size
  double CFL_nominal = enforceCFLCondition(dt_nominal);

  // create entropy viscosity
  /*
  EntropyViscosity<dim> EV(this->fe, this->n_cells, this->dof_handler, this->constraints,
      this->cell_quadrature, this->face_quadrature, this->transport_direction,
      this->cross_section_function, this->source_function, this->parameters.entropy_string,
      this->parameters.entropy_derivative_string,
      this->parameters.entropy_residual_coefficient, this->parameters.jump_coefficient,
      domain_volume, this->parameters.EV_time_discretization, low_order_viscosity,
      this->inviscid_ss_matrix, this->high_order_diffusion_matrix, this->high_order_ss_matrix);

  // create FCT object
  FCT<dim> fct(this->dof_handler, this->triangulation, lumped_mass_matrix,
      consistent_mass_matrix, this->linear_solver, this->constrained_sparsity_pattern,
      dirichlet_nodes, this->n_dofs, this->dofs_per_cell, this->parameters.do_not_limit);
      */

  // create SSP Runge-Kutta time integrator object
  SSPRKTimeIntegrator<dim> ssprk(this->parameters.time_discretization_option, this->n_dofs,
      this->linear_solver, this->constrained_sparsity_pattern);

  // time loop
  double t_new = 0.0;
  double t_old = 0.0;
  double dt = dt_nominal;
  double old_dt = dt_nominal;
  double older_dt = dt_nominal;
  old_solution = this->new_solution;
  older_solution = this->new_solution;
  oldest_solution = this->new_solution;
  const double t_end = this->parameters.end_time;
  bool in_transient = true;
  Vector<double> tmp_vector(this->n_dofs);
  unsigned int n = 1; // time step index
  while (in_transient)
  {
    // shorten time step size if new time would overshoot end time
    dt = dt_nominal;
    if (t_old + dt >= t_end)
    {
      dt = t_end - t_old;
      in_transient = false;
    }
    t_new = t_old + dt;
    std::cout << "   time step " << n << ": t = " << t_old << "->" << t_new
        << std::endl;

    // initialize SSPRK time step
    ssprk.initialize_time_step(old_solution, dt);

    switch (this->parameters.viscosity_option)
    {
      case 0 :
      { // unmodified Galerkin scheme

        takeGalerkinStep(ssprk);

        break;
      }
      case 1 :
      { // solve low-order system

        takeLowOrderStep(ssprk);

        break;
      }
      /*
      case 2 :
      { // high-order system with entropy viscosity
        // compute EV only at beginning of time step
        if (parameters.EV_time_discretization != EntropyViscosity<dim>::FE)
        {
          // recompute high-order steady-state matrix
          EV.recompute_high_order_ss_matrix(old_solution, older_solution,
              oldest_solution, old_dt, older_dt, t_old);
        }

        for (unsigned int i = 0; i < ssprk.n_stages; ++i)
        {
          // get stage time
          double t_stage = ssprk.get_stage_time();

          // recompute steady-state rhs if it is time-dependent
          if (source_time_dependent)
          {
            assemble_ss_rhs(t_stage);
          }

          if (parameters.EV_time_discretization == EntropyViscosity<dim>::FE)
          {
            // compute Galerkin solution
            ssprk.step(consistent_mass_matrix, inviscid_ss_matrix, ss_rhs,
                false);

            // get Galerkin solution
            ssprk.get_intermediate_solution(new_solution);

            // get old stage solution
            ssprk.get_stage_solution(i, old_stage_solution);

            // recompute high-order steady-state matrix
            EV.recompute_high_order_ss_matrix(new_solution, old_stage_solution,
              old_solution, dt, old_dt, t_stage);
          }

          // advance by an SSPRK step
          ssprk.step(consistent_mass_matrix, high_order_ss_matrix, ss_rhs,
              true);
        }
        // retrieve the final solution
        ssprk.get_new_solution(new_solution);

        break;
      }
      case 3 :
      { // EV FCT
        for (unsigned int i = 0; i < ssprk.n_stages; ++i)
        {
          // get stage time
          double t_stage = ssprk.get_stage_time();

          // recompute steady-state rhs if it is time-dependent
          if (source_time_dependent)
          {
            assemble_ss_rhs(t_stage);
          }

          // compute Galerkin solution
          ssprk.step(consistent_mass_matrix, inviscid_ss_matrix, ss_rhs,
              false);

          // get Galerkin solution
          ssprk.get_intermediate_solution(new_solution);

          // get old stage solution
          ssprk.get_stage_solution(i, old_stage_solution);

          // recompute high-order steady-state matrix
          EV.recompute_high_order_ss_matrix(new_solution, old_stage_solution,
              old_solution, dt, old_dt, t_stage);

          // advance by an SSPRK step
          ssprk.step(consistent_mass_matrix, high_order_ss_matrix, ss_rhs,
              false);

          // get old stage solution
          ssprk.get_stage_solution(i, old_stage_solution);

          // get intermediate solution
          ssprk.get_intermediate_solution(new_solution);

          // perform FCT
          fct.solve_FCT_system(new_solution, old_stage_solution,
              low_order_ss_matrix, ss_rhs, dt, low_order_diffusion_matrix,
              high_order_diffusion_matrix);

          // set stage solution to be FCT solution for this stage
          ssprk.set_intermediate_solution(new_solution);

          // finish computing stage solution
          ssprk.complete_stage_solution();

        }
        // retrieve the final solution
        ssprk.get_new_solution(new_solution);

        break;
      }
      case 4 :
      { // Galerkin FCT
        for (unsigned int i = 0; i < ssprk.n_stages; ++i)
        {
          // get stage time
          double t_stage = ssprk.get_stage_time();

          // recompute steady-state rhs if it is time-dependent
          if (source_time_dependent)
          {
            assemble_ss_rhs(t_stage);
          }

          // compute Galerkin solution
          ssprk.step(consistent_mass_matrix, inviscid_ss_matrix, ss_rhs,
              false);

          // get Galerkin solution
          ssprk.get_intermediate_solution(new_solution);

          // get old stage solution
          ssprk.get_stage_solution(i, old_stage_solution);

          // perform FCT
          fct.solve_FCT_system(new_solution, old_stage_solution,
              low_order_ss_matrix, ss_rhs, dt, low_order_diffusion_matrix,
              high_order_diffusion_matrix);

          // set stage solution to be FCT solution for this stage
          ssprk.set_intermediate_solution(new_solution);

          // finish computing stage solution
          ssprk.complete_stage_solution();
        }
        // retrieve the final solution
        ssprk.get_new_solution(new_solution);

        break;
      }
      */
      default :
      {
        Assert(false, ExcNotImplemented());
        break;
      }
    }

    // update old solution, time, and time step size
    oldest_solution = older_solution;
    older_solution = old_solution;
    old_solution = this->new_solution;
    older_dt = old_dt;
    old_dt = dt;
    t_old = t_new;

    // increment time step index
    n++;
  }

  // evaluate errors for convergence study
  this->postprocessor->evaluate_error(this->new_solution, this->dof_handler, *this->triangulation);

  // output grid and solution and print convergence results
  this->postprocessor->output_results(this->new_solution, this->dof_handler, *this->triangulation);
}

/**
 * Assembles the consistent and lumped mass matrices.
 */
template<int dim>
void TransientExecutioner<dim>::assembleMassMatrices()
{
  FEValues < dim
      > fe_values(this->fe, this->cell_quadrature, update_values | update_JxW_values);

  std::vector<types::global_dof_index> local_dof_indices(this->dofs_per_cell);

  FullMatrix<double> local_mass_consistent(this->dofs_per_cell, this->dofs_per_cell);
  FullMatrix<double> local_mass_lumped(this->dofs_per_cell, this->dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
      this->dof_handler.begin_active(), endc = this->dof_handler.end();
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
          local_mass_consistent(i, j) += fe_values.shape_value(i, q)
              * fe_values.shape_value(j, q) * fe_values.JxW(q);
          local_mass_lumped(i, i) += fe_values.shape_value(i, q)
              * fe_values.shape_value(j, q) * fe_values.JxW(q);
        }

    // add to global mass matrices with contraints
    this->constraints.distribute_local_to_global(local_mass_consistent,
        local_dof_indices, consistent_mass_matrix);
    this->constraints.distribute_local_to_global(local_mass_lumped,
        local_dof_indices, lumped_mass_matrix);
  }
}

/**
 * Checks CFL condition and adjusts time step size as necessary.
 *
 * @param[in,out] dt time step size for current time step
 * @return CFL number for this time step
 */
template<int dim>
double TransientExecutioner<dim>::enforceCFLCondition(double & dt)
{
  // CFL is dt*speed/dx
  double max_speed_dx = 0.0; // max(speed/dx); max over all i of A(i,i)/mL(i,i)
  for (unsigned int i = 0; i < this->n_dofs; ++i)
    max_speed_dx = std::max(max_speed_dx,
        std::abs(this->low_order_ss_matrix(i, i)) / lumped_mass_matrix(i, i));

  // compute CFL number
  double proposed_CFL = dt * max_speed_dx;

  // if computed CFL number is greater than the set limit, then adjust dt
  double adjusted_CFL = proposed_CFL;
  if (proposed_CFL > this->parameters.CFL_limit)
  {
    adjusted_CFL = this->parameters.CFL_limit;
    dt = adjusted_CFL / max_speed_dx;
  }

  // return adjusted CFL number
  return adjusted_CFL;
}

/**
 * Takes time step with the Galerkin (no viscosity) method.
 */
template<int dim>
void TransientExecutioner<dim>::takeGalerkinStep(SSPRKTimeIntegrator<dim> & ssprk)
{
  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
  {
    // get stage time
    double t_stage = ssprk.get_stage_time();

    // recompute steady-state rhs if it is time-dependent
    if (source_is_time_dependent)
    {
      this->assembleSteadyStateRHS(t_stage);
    }

    // advance by an SSPRK step
    ssprk.step(consistent_mass_matrix, this->inviscid_ss_matrix, this->ss_rhs,
        true);
  }
  // retrieve the final solution
  ssprk.get_new_solution(this->new_solution);
}

/**
 * Takes time step with the low-order viscosity method.
 */
template<int dim>
void TransientExecutioner<dim>::takeLowOrderStep(SSPRKTimeIntegrator<dim> & ssprk)
{
  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
  {
    // get stage time
    double t_stage = ssprk.get_stage_time();

    // recompute steady-state rhs if it is time-dependent
    if (source_is_time_dependent)
    {
      this->assembleSteadyStateRHS(t_stage);
    }

    // advance by an SSPRK step
    ssprk.step(lumped_mass_matrix, this->low_order_ss_matrix, this->ss_rhs,
        true);
  }
  // retrieve the final solution
  ssprk.get_new_solution(this->new_solution);
}
