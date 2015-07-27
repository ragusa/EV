/**
 * Constructor.
 */
template<int dim>
TransientExecutioner<dim>::TransientExecutioner(
  const TransportParameters<dim> & parameters,
  const Triangulation<dim> & triangulation,
  const Tensor<1, dim> & transport_direction,
  const FunctionParser<dim> & cross_section_function,
  FunctionParser<dim> & source_function, Function<dim> & incoming_function) :
    Executioner<dim>(parameters, triangulation, transport_direction,
      cross_section_function, source_function, incoming_function)
{
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
  // enforce CFL condition on nominal dt size
  CFL_nominal = enforce_CFL_condition(dt_nominal);

  // update dt displayed in convergence table
  postprocessor.update_dt(dt_nominal);

  // print dt
  std::cout << "   Nominal time step size: " << dt_nominal << std::endl;
  std::cout << "   Nominal CFL number: " << CFL_nominal << std::endl;

  // interpolate initial conditions
  initial_conditions.set_time(0.0);
  VectorTools::interpolate(dof_handler, initial_conditions, new_solution);
  constraints.distribute(new_solution);

  // if last cycle, output initial conditions if user requested
  if ((cycle == parameters.n_refinement_cycles - 1)
      and (parameters.output_initial_solution))
  {
    std::stringstream IC_filename_ss;
    IC_filename_ss << "solution_" << parameters.problem_id << "_initial";
    postprocessor.output_solution(new_solution, dof_handler,
        IC_filename_ss.str());
  }

  // create SSP Runge-Kutta time integrator object
  SSPRKTimeIntegrator<dim> ssprk(parameters.time_discretization_option, n_dofs,
      linear_solver, constrained_sparsity_pattern);

  // time loop
  double t_new = 0.0;
  double t_old = 0.0;
  double dt = dt_nominal;
  double old_dt = dt_nominal;
  double older_dt = dt_nominal;
  old_solution = new_solution;
  older_solution = new_solution;
  oldest_solution = new_solution;
  const double t_end = parameters.end_time;
  bool in_transient = true;
  Vector<double> tmp_vector(n_dofs);
  unsigned int n = 1; // time step index
  while (in_transient)
  {
    // shorten dt if new time would overshoot end time
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

    switch (parameters.viscosity_option)
    {
      case 0 :
      { // unmodified Galerkin scheme

        for (unsigned int i = 0; i < ssprk.n_stages; ++i)
        {
          // get stage time
          double t_stage = ssprk.get_stage_time();

          // recompute steady-state rhs if it is time-dependent
          if (source_time_dependent)
          {
            assemble_ss_rhs(t_stage);
          }

          // advance by an SSPRK step
          ssprk.step(consistent_mass_matrix, inviscid_ss_matrix, ss_rhs,
              true);
        }
        // retrieve the final solution
        ssprk.get_new_solution(new_solution);

        break;
      }
      case 1 :
      { // solve low-order system

        for (unsigned int i = 0; i < ssprk.n_stages; ++i)
        {
          // get stage time
          double t_stage = ssprk.get_stage_time();

          // recompute steady-state rhs if it is time-dependent
          if (source_time_dependent)
          {
            assemble_ss_rhs(t_stage);
          }

          // advance by an SSPRK step
          ssprk.step(lumped_mass_matrix, low_order_ss_matrix, ss_rhs, true);
        }
        // retrieve the final solution
        ssprk.get_new_solution(new_solution);

        break;
      }
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
            EV.recompute_high_order_ss_matrix(new_solution,
                old_stage_solution, old_solution, // not used; supply arbitrary argument
                dt, old_dt,       // not used; supply arbitrary argument
                t_stage);
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
      default :
      {
        Assert(false, ExcNotImplemented());
        break;
      }
    }

    // update old solution, time, and time step size
    oldest_solution = older_solution;
    older_solution = old_solution;
    old_solution = new_solution;
    older_dt = old_dt;
    old_dt = dt;
    t_old = t_new;

    // increment time step index
    n++;
  }
}
