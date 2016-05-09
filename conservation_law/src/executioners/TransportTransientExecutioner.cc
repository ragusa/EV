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
    source_is_time_dependent(problem_parameters_.source_is_time_dependent),
    temporal_discretization(parameters_.temporal_discretization),
    initial_conditions_function(
      &(problem_parameters_.initial_conditions_function)),
    dt_nominal(nominal_dt_)
{
  // initialize matrices
  consistent_mass_matrix.reinit(this->constrained_sparsity_pattern);
  lumped_mass_matrix.reinit(this->constrained_sparsity_pattern);

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
  ss_rhs_new.reinit(this->n_dofs);
  tmp_vector.reinit(this->n_dofs);

  // initial total iteration counts
  total_entropy_viscosity_iterations = 0;
  total_fct_iterations = 0;
}

/**
 * \brief Runs transient executioner.
 */
template <int dim>
void TransportTransientExecutioner<dim>::run()
{
  // assemble inviscid steady-state matrix
  this->assembleInviscidSteadyStateMatrix();

  // compute reaction vector
  this->compute_reaction_vector();

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

  // determine nominal time step size
  switch (this->parameters.time_step_size_option)
  {
    case TimeStepSizeOption::constant:
    {
      // leave time step as is
      break;
    }
    case TimeStepSizeOption::cfl_dmp:
    {
      // compute time step size from DMP CFL condition
      dt_nominal = compute_dt_from_dmp_cfl_condition();
      break;
    }
    default:
    {
      AssertThrow(false, ExcNotImplemented());
      break;
    }
  }

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
  unsigned int n = 0; // time step index
  bool in_transient = true;
  while (in_transient)
  {
    // increment time step index
    n++;

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

    // compute new solution in derived class
    compute_new_solution(dt, dt_old, t_old, n);

    // increment transient counter for post-processor
    this->postprocessor->increment_transient_counter();

    // update old solution, time, and time step size
    older_solution = old_solution;
    old_solution = this->new_solution;
    dt_old = dt;
    t_old = t_new;
  }

  // evaluate errors for convergence study
  this->postprocessor->evaluate_error(
    this->new_solution, this->dof_handler, *this->triangulation);

  // output grid and solution and print convergence results if in last cycle
  this->postprocessor->output_results_if_last_cycle(
    this->new_solution, this->dof_handler, *this->triangulation);

  // FCT output
  auto base_fct = get_derived_fct();
  if (base_fct)
  {
    // output final FCT bounds if specified
    if (this->parameters.output_final_fct_bounds)
      base_fct->output_bounds(*(this->postprocessor));

    // output limiter matrix if specified
    if (this->parameters.output_limiter_matrix)
      base_fct->output_limiter_matrix();
  }

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

  // evaluate average number of iterations per time step
  const double avg_ev_iterations_per_step =
    total_entropy_viscosity_iterations / n;
  const double avg_fct_iterations_per_step = total_fct_iterations / n;
  this->cout1 << std::endl;
  this->cout1 << "Total EV  iterations: " << total_entropy_viscosity_iterations
              << std::endl;
  this->cout1 << "Total FCT iterations: " << total_fct_iterations << std::endl;
  this->cout1 << "Average EV  iterations per time step: "
              << avg_ev_iterations_per_step << std::endl;
  this->cout1 << "Average FCT iterations per time step: "
              << avg_fct_iterations_per_step << std::endl;
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
 * \brief Computes time step size from the DMP CFL condition.
 *
 * \return time step size
 */
template <int dim>
double TransportTransientExecutioner<dim>::compute_dt_from_dmp_cfl_condition()
  const
{
  // CFL is dt*speed/dx
  double max_speed_dx = 0.0; // max(speed/dx); max over all i of A(i,i)/mL(i,i)
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    max_speed_dx = std::max(max_speed_dx,
                            std::abs(this->low_order_ss_matrix(i, i)) /
                              lumped_mass_matrix(i, i));
  }

  // compute dt from DMP CFL condition
  const double dt = this->parameters.cfl / max_speed_dx;

  return dt;
}
