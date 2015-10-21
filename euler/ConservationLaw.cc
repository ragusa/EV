/**
 * \file ConservationLaw.cc
 * \brief Provides the function definitions for the ConservationLaw class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] params conservation law parameters
 */
template <int dim>
ConservationLaw<dim>::ConservationLaw(
  const ConservationLawParameters<dim> & params)
  : parameters(params),
    n_components(params.n_components),
    mapping(),
    fe(FE_Q<dim>(params.degree), params.n_components),
    dofs_per_cell(fe.dofs_per_cell),
    faces_per_cell(GeometryInfo<dim>::faces_per_cell),
    dof_handler(triangulation),
    n_q_points_per_dim(params.n_quadrature_points),
    cell_quadrature(n_q_points_per_dim),
    face_quadrature(n_q_points_per_dim),
    n_q_points_cell(cell_quadrature.size()),
    n_q_points_face(face_quadrature.size()),
    initial_conditions_strings(params.n_components),
    initial_conditions_function(params.n_components),
    exact_solution_strings(params.n_components),
    exact_solution_function()
{
}

/**
 * \brief Destructor.
 */
template <int dim>
ConservationLaw<dim>::~ConservationLaw()
{
}

/**
 * \brief Runs the problem.
 */
template <int dim>
void ConservationLaw<dim>::run()
{
  // initialize system
  initialize_system();

  // create post-processor object
  PostProcessor<dim> postprocessor(parameters,
                                   end_time,
                                   has_exact_solution,
                                   exact_solution_function,
                                   problem_name,
                                   component_names,
                                   component_interpretations,
                                   triangulation);

  // loop over adaptive refinement cycles
  for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; ++cycle)
  {
    // set cycle for post-processor
    if (cycle == parameters.n_refinement_cycles - 1)
      postprocessor.set_cycle(cycle);

    // refine mesh if not the first cycle
    if (cycle > 0)
      refine_mesh();

    // setup system; to be applied after each refinement
    setup_system();

    // print information
    std::cout << std::endl
              << "Cycle " << cycle << " of " << parameters.n_refinement_cycles - 1
              << ":" << std::endl;
    std::cout << "  Number of active cells: " << triangulation.n_active_cells()
              << std::endl;

    // interpolate the initial conditions to the grid, and apply constraints
    VectorTools::interpolate(
      dof_handler, initial_conditions_function, new_solution);
    apply_Dirichlet_BC(0.0);
    constraints.distribute(new_solution);

    // set old solution to the current solution
    old_solution = new_solution;

    // output initial solution
    postprocessor.output_solution(new_solution, dof_handler, "solution_initial");
    postprocessor.output_solution_transient(
      new_solution, dof_handler, "solution");

    // solve transient with selected time integrator
    switch (parameters.temporal_integrator)
    {
      case ConservationLawParameters<dim>::runge_kutta: // Runge-Kutta
        solve_runge_kutta(postprocessor);
        break;
      default:
        Assert(false, ExcNotImplemented());
        break;
    }

    // evaluate errors for convergence study
    postprocessor.evaluate_error(new_solution, dof_handler, triangulation);

  } // end of refinement loop

  // output grid and solution and print convergence results
  output_results(postprocessor);

  // output viscosity if requested
  if (parameters.output_viscosity)
  {
    // output final viscosities if non-constant viscosity used
    switch (parameters.viscosity_type)
    {
      case ConservationLawParameters<dim>::none:
        break;
      case ConservationLawParameters<dim>::constant:
        break;
      case ConservationLawParameters<dim>::old_first_order:
        postprocessor.output_cell_map(
          first_order_viscosity, "loworderviscosity", dof_handler);
        break;
      case ConservationLawParameters<dim>::max_principle:
        postprocessor.output_cell_map(
          first_order_viscosity, "loworderviscosity", dof_handler);
        break;
      case ConservationLawParameters<dim>::entropy:
        postprocessor.output_cell_map(
          first_order_viscosity, "loworderviscosity", dof_handler);

        // if low-order viscosity is used for first time step and only 1 time
        // step is taken, then the size of the entropy viscosity map is 0,
        // so that map cannot be output
        if (entropy_viscosity.size() > 0)
          postprocessor.output_cell_map(
            entropy_viscosity, "entropyviscosity", dof_handler);

        postprocessor.output_cell_map(
          viscosity, "highorderviscosity", dof_handler);
        break;
      default:
        Assert(false, ExcNotImplemented());
        break;
    }
  }
}

/**
 * \brief Initializes system.
 */
template <int dim>
void ConservationLaw<dim>::initialize_system()
{
  // get component names and interpretations
  component_names = get_component_names();
  component_interpretations = get_component_interpretations();

  // define problem parameters
  define_problem();

  // make initial triangulation
  triangulation.refine_global(parameters.initial_refinement_level);

  // determine end time if chose to use default end time
  if (parameters.use_default_end_time && has_default_end_time)
    end_time = default_end_time;
  else
    end_time = parameters.end_time;

  // create constants used for parsed functions
  constants["pi"] = numbers::PI;

  // create string of variables to be used in function parser objects
  std::string variables = FunctionParser<dim>::default_variable_names();
  /*
  std::string variables;
  if (dim == 1)
    variables = "x";
  else if (dim == 2)
    variables = "x,y";
  else if (dim == 3)
    variables = "x,y,z";
  else
    Assert(false, ExcInvalidState());
  */

  // add t for time-dependent function parser objects
  std::string time_dependent_variables = variables + ",t";

  // initialize Dirichlet boundary functions if needed
  if (boundary_conditions_type == "dirichlet" &&
      !(use_exact_solution_as_dirichlet_bc))
  {
    dirichlet_function.resize(n_boundaries);
    for (unsigned int boundary = 0; boundary < n_boundaries; ++boundary)
    {
      dirichlet_function[boundary] = new FunctionParser<dim>(n_components);
      dirichlet_function[boundary]->initialize(
        time_dependent_variables,
        dirichlet_function_strings[boundary],
        constants,
        true);
    }
  }

  // initialize initial conditions function
  initial_conditions_function.initialize(
    variables, initial_conditions_strings, constants, false);

  // initialize Runge-Kutta data if RK is being used
  if (parameters.temporal_integrator ==
      ConservationLawParameters<dim>::runge_kutta)
    initialize_runge_kutta();
}

/** \brief Assigns Butcher tableau constants.
 */
template <int dim>
void ConservationLaw<dim>::initialize_runge_kutta()
{
  // get RK parameters a, b, and c (Butcher tableau)
  rk.s = 0;
  switch (parameters.time_discretization)
  {
    case ConservationLawParameters<dim>::ERK1:
      rk.s = 1;
      break;
    case ConservationLawParameters<dim>::ERK2:
      rk.s = 2;
      break;
    case ConservationLawParameters<dim>::ERK3:
      rk.s = 3;
      break;
    case ConservationLawParameters<dim>::ERK4:
      rk.s = 4;
      break;
    case ConservationLawParameters<dim>::SDIRK22:
      rk.s = 2;
      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // allocate memory for constants
  rk.a.resize(rk.s);
  for (int i = 0; i < rk.s; ++i)
    rk.a[i].reinit(rk.s);
  rk.b.resize(rk.s);
  rk.c.resize(rk.s);

  // assign constants
  double gamma, sigma;
  switch (parameters.time_discretization)
  {
    case ConservationLawParameters<dim>::ERK1:
      rk.b[0] = 1;
      rk.c[0] = 0;
      rk.is_explicit = true;
      break;
    case ConservationLawParameters<dim>::ERK2:
      rk.a[1][0] = 0.5;
      rk.b[0] = 0;
      rk.b[1] = 1;
      rk.c[0] = 0;
      rk.c[1] = 0.5;
      rk.is_explicit = true;
      break;
    case ConservationLawParameters<dim>::ERK3:
      rk.a[1][0] = 1.0;
      rk.a[2][0] = 0.25;
      rk.a[2][1] = 0.25;
      rk.b[0] = 1. / 6;
      rk.b[1] = 1. / 6;
      rk.b[2] = 4. / 6;
      rk.c[0] = 0;
      rk.c[1] = 1.0;
      rk.c[2] = 0.5;
      rk.is_explicit = true;
      break;
    case ConservationLawParameters<dim>::ERK4:
      rk.a[1][0] = 0.5;
      rk.a[2][0] = 0;
      rk.a[2][1] = 0.5;
      rk.a[3][0] = 0;
      rk.a[3][1] = 0;
      rk.a[3][2] = 1;
      rk.b[0] = 1. / 6;
      rk.b[1] = 1. / 3;
      rk.b[2] = 1. / 3;
      rk.b[3] = 1. / 6;
      rk.c[0] = 0;
      rk.c[1] = 0.5;
      rk.c[2] = 0.5;
      rk.c[3] = 1;
      rk.is_explicit = true;
      break;
    case ConservationLawParameters<dim>::SDIRK22:
      gamma = 1.0 - 1.0 / std::sqrt(2.0);
      sigma = 1.0 - gamma;
      rk.a[0][0] = gamma;
      rk.a[1][0] = sigma;
      rk.a[1][1] = gamma;
      rk.b[0] = sigma;
      rk.b[1] = gamma;
      rk.c[0] = gamma;
      rk.c[1] = 1.0;
      rk.is_explicit = false;
      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // allocate stage vectors for steady-state residual evaluations
  rk.f.resize(rk.s);

  // test to see if last row of A matrix is the same as the b vector;
  // this implies the last stage solution is the new solution and thus
  // the final linear combination is not required.
  rk.solution_computed_in_last_stage = true;
  for (int i = 0; i < rk.s; ++i)
    if (std::abs(rk.a[rk.s - 1][i] - rk.b[i]) > 1.0e-30)
      rk.solution_computed_in_last_stage = false;
}

/** \brief Computes error for adaptive mesh refinement for a time
 *         step and adds it to an error sum for all time steps.
 */
template <int dim>
void ConservationLaw<dim>::compute_error_for_refinement()
{
  Vector<float> estimated_error_per_cell_time_step(
    triangulation.n_active_cells());

  KellyErrorEstimator<dim>::estimate(
    dof_handler,
    face_quadrature,
    // for now, assume no Neumann boundary conditions,
    //  so the following argument may be empty
    typename FunctionMap<dim>::type(),
    new_solution,
    estimated_error_per_cell_time_step);

  estimated_error_per_cell += estimated_error_per_cell_time_step;
}

/** \brief Adaptively refines mesh based on estimated error per cell.
 */
template <int dim>
void ConservationLaw<dim>::refine_mesh()
{
  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                  estimated_error_per_cell,
                                                  parameters.refinement_fraction,
                                                  parameters.coarsening_fraction);

  triangulation.execute_coarsening_and_refinement();
}

/** \brief Sets up the system before solving.
 *
 *         This function is to be applied after each refinement. It
 *         allocates memory, sets up constraints, makes the sparsity pattern,
 *         and reinitializes the system matrix with the sparsity pattern.
 */
template <int dim>
void ConservationLaw<dim>::setup_system()
{
  // clear maps
  cell_diameter.clear();
  max_flux_speed_cell.clear();
  viscosity.clear();
  first_order_viscosity.clear();
  entropy_viscosity.clear();

  // clear and distribute dofs
  dof_handler.clear();
  dof_handler.distribute_dofs(fe);
  n_dofs = dof_handler.n_dofs();
  std::cout << "Number of degrees of freedom: " << n_dofs << std::endl;

  // determine Dirichlet nodes
  get_dirichlet_nodes();

  // update cell sizes and minimum cell size
  update_cell_sizes();

  // make constraints
  constraints.clear();
  // hanging node constraints
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  // Dirichlet contraints (for t = 0)
  if (boundary_conditions_type == "dirichlet")
  {
    for (unsigned int boundary = 0; boundary < n_boundaries; ++boundary)
      for (unsigned int component = 0; component < n_components; ++component)
      {
        // specify to impose Dirichlet BC for all components of solution
        std::vector<bool> component_mask(n_components, true);

        if (use_exact_solution_as_dirichlet_bc)
        {
          exact_solution_function->set_time(0.0);
          VectorTools::interpolate_boundary_values(dof_handler,
                                                   boundary,
                                                   *exact_solution_function,
                                                   constraints,
                                                   component_mask);
        }
        else
        {
          dirichlet_function[boundary]->set_time(0.0);
          VectorTools::interpolate_boundary_values(
            dof_handler,
            boundary,
            *(dirichlet_function[boundary]),
            constraints,
            component_mask);
        }
      }
  }
  constraints.close();

  // create sparsity patterns
  DynamicSparsityPattern dynamic_constrained_sparsity_pattern(n_dofs);
  DynamicSparsityPattern dynamic_unconstrained_sparsity_pattern(n_dofs);
  DoFTools::make_sparsity_pattern(
    dof_handler, dynamic_constrained_sparsity_pattern, constraints, false);
  DoFTools::make_sparsity_pattern(dof_handler,
                                  dynamic_unconstrained_sparsity_pattern);
  constrained_sparsity_pattern.copy_from(dynamic_constrained_sparsity_pattern);
  unconstrained_sparsity_pattern.copy_from(
    dynamic_unconstrained_sparsity_pattern);

  // reinitialize matrices with sparsity pattern
  consistent_mass_matrix.reinit(constrained_sparsity_pattern);
  lumped_mass_matrix.reinit(constrained_sparsity_pattern);
  system_matrix.reinit(constrained_sparsity_pattern);

  // if using maximum-principle preserving definition of first-order viscosity,
  // then compute bilinear forms and viscous fluxes
  if (parameters.viscosity_type == ConservationLawParameters<dim>::max_principle)
  {
    viscous_bilinear_forms.reinit(unconstrained_sparsity_pattern);
    compute_viscous_bilinear_forms();
    viscous_fluxes.reinit(unconstrained_sparsity_pattern);
  }

  // assemble mass matrix
  assemble_mass_matrix();

  // resize vectors
  old_solution.reinit(n_dofs);
  new_solution.reinit(n_dofs);
  exact_solution.reinit(n_dofs);
  solution_step.reinit(n_dofs);
  system_rhs.reinit(n_dofs);
  max_values.reinit(n_dofs);
  min_values.reinit(n_dofs);
  estimated_error_per_cell.reinit(triangulation.n_active_cells());

  // allocate memory for steady-state residual evaluations if Runge-Kutta is
  // used
  if (parameters.temporal_integrator ==
      ConservationLawParameters<dim>::runge_kutta)
    for (int i = 0; i < rk.s; ++i)
      rk.f[i].reinit(n_dofs);

  // perform additional setup required by derived classes
  perform_additional_setup();
}

/** \brief Updates the cell sizes map and minimum cell size.
 */
template <int dim>
void ConservationLaw<dim>::update_cell_sizes()
{
  // fill cell size map and find minimum cell size
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  // reset minimum cell size to an arbitrary cell size such as the first cell
  minimum_cell_diameter = cell->diameter();

  // update the cell diameters and minimum cell diameter
  for (; cell != endc; ++cell)
  {
    cell_diameter[cell] = cell->diameter();
    minimum_cell_diameter = std::min(minimum_cell_diameter, cell_diameter[cell]);
  }
}

/** \brief Assembles the mass matrices and applies constraints.
 */
template <int dim>
void ConservationLaw<dim>::assemble_mass_matrix()
{
  // assemble lumped mass matrix
  //----------------------------
  assemble_lumped_mass_matrix();

  // assemble consistent mass matrix
  //----------------------------
  Function<dim> * dummy_function = 0;
  MatrixTools::create_mass_matrix(dof_handler,
                                  cell_quadrature,
                                  consistent_mass_matrix,
                                  dummy_function,
                                  constraints);

  // output mass matrix
  if (parameters.output_mass_matrix)
  {
    // output lumped mass matrix
    std::ofstream lumped_mass_matrix_out("output/lumped_mass_matrix.txt");
    lumped_mass_matrix.print_formatted(
      lumped_mass_matrix_out, 10, true, 0, "0", 1);
    lumped_mass_matrix_out.close();

    // output consistent mass matrix
    std::ofstream consistent_mass_matrix_out("output/consistent_mass_matrix.txt");
    consistent_mass_matrix.print_formatted(
      consistent_mass_matrix_out, 10, true, 0, "0", 1);
    consistent_mass_matrix_out.close();
  }
}

/** \brief Applies Dirichlet boundary conditions at a particular time.
 *
 *  \param[in] time current time; Dirichlet BC may be time-dependent.
 */
template <int dim>
void ConservationLaw<dim>::apply_Dirichlet_BC(const double & time)
{
  // map of global dof ID to boundary value, to be computed using provided
  // function
  std::map<unsigned int, double> boundary_values;

  if (boundary_conditions_type == "dirichlet")
  {
    // loop over boundary IDs
    for (unsigned int boundary = 0; boundary < n_boundaries; ++boundary)
      // loop over components
      for (unsigned int component = 0; component < n_components; ++component)
      {
        // mask other components
        std::vector<bool> component_mask(n_components, false);
        component_mask[component] = true;

        // fill boundary_values with boundary values
        if (use_exact_solution_as_dirichlet_bc)
        {
          exact_solution_function->set_time(time);
          VectorTools::interpolate_boundary_values(dof_handler,
                                                   boundary,
                                                   *exact_solution_function,
                                                   boundary_values,
                                                   component_mask);
        }
        else
        {
          dirichlet_function[boundary]->set_time(time);
          VectorTools::interpolate_boundary_values(
            dof_handler,
            boundary,
            *(dirichlet_function[boundary]),
            boundary_values,
            component_mask);
        }
      }
  }
  // apply boundary values to the solution
  for (std::map<unsigned int, double>::const_iterator it =
         boundary_values.begin();
       it != boundary_values.end();
       ++it)
    new_solution(it->first) = (it->second);
}

/** \brief Outputs a mapped quantity at all quadrature points.
 *  \param map map between cell iterator and vector of values at quadrature
 *         points within cell.
 *  \param output_filename_base string which forms the base (without extension)
 *         of the output file.
 */
/*
template <int dim>
void ConservationLaw<dim>::output_map(
  const cell_vector_map & map, const std::string & output_filename_base) const
{
  if (dim == 1)
  {
    // get data from maps into vector of point-data pairs
    unsigned int n_cells = triangulation.n_active_cells();
    unsigned int total_n_q_points = n_cells * n_q_points_cell;
    std::vector<std::pair<double, double>> profile(total_n_q_points);

    FEValues<dim> fe_values(fe, cell_quadrature, update_quadrature_points);
    unsigned int i = 0;
    cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      for (unsigned int q = 0; q < n_q_points_cell; ++q)
      {
        const Point<dim> q_point = fe_values.quadrature_point(q);
        profile[i] = std::make_pair(q_point(0), map[cell](q));
        ++i;
      }
    }

    // sort data by quadrature point
    std::sort(profile.begin(), profile.end());

    // output vector to file
    std::ofstream output;
    std::string output_file = "output/" + output_filename_base + ".csv";
    output.open(output_file.c_str(), std::ios::out);
    for (i = 0; i < total_n_q_points; ++i)
      output << profile[i].first << "," << profile[i].second << std::endl;
    output.close();
  }
}
*/

/** \brief Outputs a mapped quantity at all cells, not all quadrature points
 * within cells
 *  \param map map between cell iterator and vector of values at quadrature
 * points
 * within cell.
 *  \param output_filename_base string which forms the base (without extension)
 * of
 * the output file.
 */
/*
template <int dim>
void ConservationLaw<dim>::output_map(
  const cell_double_map & map, const std::string & output_filename_base) const
{
  if (dim == 1)
  {
    // get data from maps into vector of point-data pairs
    unsigned int n_cells = triangulation.n_active_cells();
    unsigned int total_n_q_points = n_cells * n_q_points_cell;
    std::vector<std::pair<double, double>> profile(total_n_q_points);

    FEValues<dim> fe_values(fe, cell_quadrature, update_quadrature_points);
    unsigned int i = 0;
    cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      for (unsigned int q = 0; q < n_q_points_cell; ++q)
      {
        const Point<dim> q_point = fe_values.quadrature_point(q);
        profile[i] = std::make_pair(q_point(0), map[cell]);
        ++i;
      }
    }

    // sort data by quadrature point
    std::sort(profile.begin(), profile.end());

    // output vector to file
    std::ofstream output;
    std::string output_file = "output/" + output_filename_base + ".csv";
    output.open(output_file.c_str(), std::ios::out);
    for (i = 0; i < total_n_q_points; ++i)
      output << profile[i].first << "," << profile[i].second << std::endl;
    output.close();
  }
}
*/

/**
 * \brief Solves transient using a Runge-Kutta scheme.
 *
 * This function contains the transient loop and solves the
 * transient using explicit Runge-Kutta:
 * \f[
 *   \mathbf{M} \mathbf{y}^{n+1} = \mathbf{M} \mathbf{y}^n
 *   + \Delta t\sum\limits^s_{i=1}b_i \mathbf{r}_i
 * \f]
 * where
 * \f[
 *   \mathbf{r}_i = \mathbf{r}(t^n + c_i \Delta t, \mathbf{Y}_i)
 * \f]
 * and \f$\mathbf{Y}_i\f$ is computed from the linear solve
 * \f[
 *   \mathbf{M} \mathbf{Y}_i = \mathbf{M} \mathbf{y}^n
 *   + \Delta t\sum\limits^{i-1}_{j=1}a_{i,j} \mathbf{r}_i
 * \f]
 *
 * \param[inout] postprocessor Post-processor for outputting transient
 */
template <int dim>
void ConservationLaw<dim>::solve_runge_kutta(PostProcessor<dim> & postprocessor)
{
  old_time = 0.0;
  unsigned int n = 1; // time step index
  double t_end = end_time;
  bool in_transient = true;
  bool DMP_satisfied = true;
  while (in_transient)
  {
    // update max speed for use in CFL computation
    update_flux_speeds();

    // compute dt
    double dt;
    switch (parameters.time_step_size_method)
    {
      case ConservationLawParameters<dim>::constant_dt:
        dt = parameters.time_step_size;
        break;
      case ConservationLawParameters<dim>::cfl_condition:
        dt = compute_dt_from_cfl_condition();
        break;
      default:
        Assert(false, ExcNotImplemented());
        break;
    }
    // check end of transient and shorten last time step if necessary
    if ((old_time + dt) >= t_end)
    {
      dt = t_end - old_time;
      in_transient = false;
    }

    std::cout << "  time step " << n << ": t = "
              << "\x1b[1;34m" << old_time + dt << "\x1b[0m" << std::endl;

    update_viscosities(dt, n);

    // update old_solution to new_solution for next time step;
    // this is not done at the end of the previous time step because
    // of the time derivative term in update_viscosities() above
    old_solution = new_solution;

    // solve here
    // compute each f_i
    for (int i = 0; i < rk.s; ++i)
    {
      std::cout << "    stage " << i + 1 << " of " << rk.s << std::endl;
      // compute stage time
      current_time = old_time + rk.c[i] * dt;

      if (rk.is_explicit) // explicit Runge-Kutta
      {
        /* compute RHS of
         *
         * M*Y_i = M*y^n + sum_{j=1}^{i-1} dt*a_{i,j}*f_j
         *
         */
        system_rhs = 0.0;
        lumped_mass_matrix.vmult(system_rhs, old_solution);
        for (int j = 0; j < i; ++j)
          system_rhs.add(dt * rk.a[i][j], rk.f[j]);

        // solve system M*Y_i = RHS
        linear_solve(lumped_mass_matrix, system_rhs, new_solution);

        // ordinarily, Dirichlet BC need not be reapplied, but in general,
        // the Dirichlet BC can be time-dependent
        apply_Dirichlet_BC(current_time);

        // distribute hanging node constraints
        constraints.distribute(new_solution);
      }
      else
      { // implicit Runge-Kutta

        Assert(false, ExcNotImplemented());
        /*
                    // Newton solve
                    new_solution = old_solution;
                    // compute initial negative of transient residual: -F(y)
                    compute_tr_residual(i,dt);
                    // compute initial norm of transient residual - used in
                    // convergence tolerance
                    double residual_norm = system_rhs.l2_norm();
                    // compute nonlinear tolerance based on provided ATOL and
           RTOL
                    double nonlinear_tolerance = parameters.nonlinear_atol +
                       parameters.nonlinear_rtol * residual_norm;
                    // initialize convergence flag
                    bool converged = false;
                    // begin Newton loop
                    for (unsigned int iteration = 0;
                      iteration < parameters.max_nonlinear_iterations;
           ++iteration)
                    {
                       // compute steady-state Jacobian and store in
           system_matrix
                       compute_ss_jacobian();
                       // compute transient Jacobian and store in system_matrix
                       system_matrix *= rk.a[i][i]*dt;
                       system_matrix.add(-1.0,lumped_mass_matrix);
                       //compute_tr_Jacobian();
                       // Solve for Newton step
                       linear_solve(system_matrix, system_rhs, solution_step);
                       // update solution
                       new_solution += solution_step;
                       // ordinarily, Dirichlet BC need not be reapplied, but in
           this case,
                       // the Dirichlet BC can be time-dependent
                       apply_Dirichlet_BC(current_time);
                       // compute negative of transient residual: -F(y)
                       compute_tr_residual(i,dt);
                       // compute norm of transient residual
                       residual_norm = system_rhs.l2_norm();
                       std::cout << "         nonlinear iteration " << iteration
                         << ": residual norm = " << residual_norm << std::endl;
                       // check convergence
                       if (residual_norm < nonlinear_tolerance)
                       {
                          // set convergence flag
                          converged = true;
                          // break out of Newton loop
                          break;
                       }
                    }
                    // exit program if solution did not converge
                    if (not converged) {
                       std::cout << "Solution did not converge within maximum
           number of
                         nonlinear iterations."
                          << " Program terminated." << std::endl;
                       std::exit(1);
                    }
        */
      }

      // residual from solution of previous step is reused in first stage
      // (unless this is the first time step)
      if ((n == 1) || (i != 0))
      {
        compute_ss_residual(rk.f[i]);
      }
    }

    /* compute the solution using the computed f_i's:
     *
     * M*y^{n+1} = M*y^n + dt*sum_{i=1}^s b_i*f_i
     */
    // if the next time step solution was already computed in the last stage,
    // then there is no need to perform the linear combination of f's. Plus,
    // since f[0] is always computed from the previous time step solution,
    // re-use the end f as f[0] for the next time step
    if (rk.solution_computed_in_last_stage)
    {
      // end of step residual can be used as beginning of step residual for next
      // step
      rk.f[0] = rk.f[rk.s - 1];
    }
    else
    {
      // solve M*y^{n+1} = M*y^n + dt*sum_{i=1}^s b_i*f_i
      system_rhs = 0.0;
      lumped_mass_matrix.vmult(system_rhs, old_solution);
      for (int i = 0; i < rk.s; ++i)
        system_rhs.add(dt * rk.b[i], rk.f[i]);
      linear_solve(lumped_mass_matrix, system_rhs, new_solution);

      // apply Dirichlet constraints
      apply_Dirichlet_BC(current_time);

      // distribute hanging node constraints
      constraints.distribute(new_solution);

      // end of step residual can be used as beginning of step residual for next
      // step
      compute_ss_residual(rk.f[0]);
    }

    // output solution transient if specified
    postprocessor.output_solution_transient(
      new_solution, dof_handler, "solution");

    // compute max principle min and max values
    compute_max_principle_quantities();

    // increment time
    old_time += dt;
    current_time = old_time; // used in exact solution function in output_solution
    n++;

    // check that there are no NaNs in solution
    check_nan();

    // check for steady-state
    Vector<double> solution_change = new_solution;
    solution_change.add(-1.0, old_solution);
    double solution_change_norm = solution_change.l2_norm();
    double solution_normalization = old_solution.l2_norm();
    double steady_state_error;
    if (std::abs(solution_normalization) < 1.0e-15)
      steady_state_error = solution_change_norm;
    else
      steady_state_error = solution_change_norm / solution_normalization;
    std::cout << "    Steady-state error = "
              << "\x1b[1;35m" << steady_state_error << "\x1b[0m" << std::endl;
    if (steady_state_error < parameters.steady_state_tolerance)
    {
      in_transient = false;
      std::cout << "\x1b[1;32m"
                << "  Converged to steady-state."
                << "\x1b[0m" << std::endl;
    }

    // check that DMP is satisfied at all time steps
    if (parameters.viscosity_type ==
        ConservationLawParameters<dim>::max_principle)
    {
      bool DMP_satisfied_this_time_step = check_DMP(n);
      DMP_satisfied = DMP_satisfied and DMP_satisfied_this_time_step;
    }

    // compute error for adaptive mesh refinement
    compute_error_for_refinement();

  } // end of time loop

  // report if DMP was satisfied at all time steps
  if (parameters.viscosity_type == ConservationLawParameters<dim>::max_principle)
  {
    if (DMP_satisfied)
      std::cout << "DMP was satisfied at all time steps" << std::endl;
    else
      std::cout << "DMP was NOT satisfied at all time steps" << std::endl;
  }
}

/** \brief Computes time step size using the CFL condition
 *
 *  The CFL condition for stability is the following:
 *  \f[
 *    \nu = \frac{\lambda_{max}\Delta t}{h_{min}}\le 1 ,
 *  \f]
 *  where \f$\lambda_{max}\f$ is the maximum speed in the domain,
 *  \f$\Delta t\f$ is the time step size, and \f$h_{min}\f$
 *  is the minimum mesh size in the domain. The user supplies the
 *  CFL number \f$\nu^{user}\f$ (which should be less than 1),
 *  and the time step is calculated as
 *  \f[
 *    \Delta t = \frac{\nu^{user}h_{min}}{\lambda_{max}} .
 *  \f]
 */
template <int dim>
double ConservationLaw<dim>::compute_dt_from_cfl_condition()
{
  return parameters.cfl * minimum_cell_diameter / max_flux_speed;
}

/** \brief Computes the CFL number.
 *
 *  The CFL number is the following:
 *  \f[
 *    \nu = \frac{\lambda_{max}\Delta t}{h_{min}} ,
 *  \f]
 *  where \f$\lambda_{max}\f$ is the maximum speed in the domain,
 *  \f$\Delta t\f$ is the time step size, and \f$h_{min}\f$
 *  is the minimum mesh size in the domain.
 *
 *  \param[in] dt time step size
 *  \return CFL number
 */
template <int dim>
double ConservationLaw<dim>::compute_cfl_number(const double & dt) const
{
  return dt * max_flux_speed / minimum_cell_diameter;
}

/** \brief Adds the viscous bilinear form for maximum-principle preserving
 * viscosity
 */
template <int dim>
void ConservationLaw<dim>::add_maximum_principle_viscosity_bilinear_form(
  Vector<double> & f)
{
  FEValues<dim> fe_values(
    fe, cell_quadrature, update_values | update_gradients | update_JxW_values);

  // allocate memory needed for cell residual and aggregation into global
  // residual
  Vector<double> cell_residual(dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reset cell residual
    cell_residual = 0;

    // add viscous bilinear form
    double cell_volume = cell->measure();
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      // compute b_K(u,\varphi_i)
      double b_i = 0.0;
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        double b_K; // local viscous bilinear form
        if (j == i)
          b_K = cell_volume;
        else
          b_K = -1.0 / (dofs_per_cell - 1.0) * cell_volume;

        b_i += new_solution(local_dof_indices[j]) * b_K;
      }
      // add viscous term for dof i; note 0 is used for quadrature point
      // because viscosity is the same for all quadrature points on a cell
      cell_residual(i) += viscosity[cell] * b_i;
    }

    // aggregate local residual into global residual
    cell->get_dof_indices(local_dof_indices);
    constraints.distribute_local_to_global(cell_residual, local_dof_indices, f);
  } // end cell loop
}

/**
 * \brief Solves the linear system \f$A x = b\f$.
 *
 * \param[in] A the system matrix
 * \param[in] b the right-hand side vector
 * \param[out] x the solution vector
 */
template <int dim>
void ConservationLaw<dim>::linear_solve(const SparseMatrix<double> & A,
                                        const Vector<double> & b,
                                        Vector<double> & x)
{
  switch (parameters.linear_solver)
  {
    // UMFPACK
    case ConservationLawParameters<dim>::direct:
    {
      SparseDirectUMFPACK A_umfpack;
      A_umfpack.initialize(A);
      A_umfpack.vmult(x, b);
      break;
    }
    // GMRes
    case ConservationLawParameters<dim>::gmres:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
    // CG with SSOR
    case ConservationLawParameters<dim>::cg:
    {
      SolverControl solver_control(parameters.max_linear_iterations,
                                   parameters.linear_atol);
      SolverCG<> solver(solver_control);

      PreconditionSSOR<> preconditioner;
      // initialize. The second parameter is a relaxation parameter
      // between 1 and 2.
      preconditioner.initialize(A, 1.2);

      solver.solve(A, x, b, preconditioner);
      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }

  // distribute constraints
  constraints.distribute(x);
}

/**
 * \brief Updates viscosity at each quadrature point in each cell.
 *
 * \param[in] dt time step size
 * \param[in] n time index
 */
template <int dim>
void ConservationLaw<dim>::update_viscosities(const double & dt,
                                              const unsigned int & n)
{
  switch (parameters.viscosity_type)
  {
    // no viscosity
    case ConservationLawParameters<dim>::none:
    {
      cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
        viscosity[cell] = 0.0;
      break;
    }
    // constant viscosity
    case ConservationLawParameters<dim>::constant:
    {
      cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
        viscosity[cell] = parameters.constant_viscosity_value;
      break;
    }
    // old first order viscosity
    case ConservationLawParameters<dim>::old_first_order:
    {
      update_old_low_order_viscosity(true);
      viscosity = first_order_viscosity;
      break;
    }
    // max principle viscosity
    case ConservationLawParameters<dim>::max_principle:
    {
      update_max_principle_viscosity();
      viscosity = first_order_viscosity;
      break;
    }
    // entropy viscosity
    case ConservationLawParameters<dim>::entropy:
    {
      cell_iterator cell, endc = dof_handler.end();
      if (parameters.use_low_order_viscosity_for_first_time_step && n == 1)
      {
        // compute low-order viscosities
        update_old_low_order_viscosity(true);

        // use low-order viscosities for first time step
        for (cell = dof_handler.begin_active(); cell != endc; ++cell)
          viscosity[cell] = first_order_viscosity[cell];
      }
      else
      {
        // compute low-order viscosities
        update_old_low_order_viscosity(false);

        // compute high-order viscosities
        update_entropy_viscosities(dt);
        for (cell = dof_handler.begin_active(); cell != endc; ++cell)
          viscosity[cell] =
            std::min(first_order_viscosity[cell], entropy_viscosity[cell]);
      }
      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }
}

/**
 * \brief Computes low-order viscosity for each cell.
 *
 * The low-order viscosity is computed as
 * \f[
 *   \nu_K^L = c_{max} h_K \lambda_{K,max} \,,
 * \f]
 * where \f$\lambda_{K,max}\f$ is the maximum flux speed on cell \f$K\f$.
 *
 * \param[in] use_low_order_scheme flag that the low-order viscosities are
 *            are used in a low-order scheme as opposed to an entropy
 *            viscosity scheme
 */
template <int dim>
void ConservationLaw<dim>::update_old_low_order_viscosity(const bool &)
{
  double c_max = parameters.first_order_viscosity_coef;

  // loop over cells to compute first order viscosity at each quadrature point
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    first_order_viscosity[cell] =
      std::abs(c_max * cell_diameter[cell] * max_flux_speed_cell[cell]);
  }
}

/**
 * \brief Computes the maximum-principle preserving first order viscosity
 *        for each cell.
 */
template <int dim>
void ConservationLaw<dim>::update_max_principle_viscosity()
{
  // compute viscous fluxes
  compute_viscous_fluxes();

  // local dof indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells to compute first order viscosity at each quadrature point
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    first_order_viscosity[cell] = 0.0;
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        if (i != j)
        {
          first_order_viscosity[cell] = std::max(
            first_order_viscosity[cell],
            std::abs(viscous_fluxes(local_dof_indices[i], local_dof_indices[j])) /
              (-viscous_bilinear_forms(local_dof_indices[i],
                                       local_dof_indices[j])));
        }
      }
    }
  }
}

/** \brief Computes viscous fluxes, to be used in the computation of
 *         maximum-principle preserving first order viscosity.
 *
 *         Each element of the resulting matrix, \f$V_{i,j}\f$ is computed as
 *         follows:
 *         \f[
 *            V_{i,j} =
 *\int_{S_{ij}}(\mathbf{f}'(u)\cdot\nabla\varphi_j)\varphi_i
 * d\mathbf{x}
 *         \f]
 */
template <int dim>
void ConservationLaw<dim>::compute_viscous_fluxes()
{
  viscous_fluxes = 0; // zero out matrix

  // local dof indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  FEValues<dim> fe_values(
    fe, cell_quadrature, update_values | update_gradients | update_JxW_values);
  std::vector<double> solution_values(n_q_points_cell);
  std::vector<Tensor<1, dim>> dfdu(n_q_points_cell);

  // loop over cells
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    fe_values.reinit(cell);
    fe_values.get_function_values(new_solution, solution_values);

    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // add local viscous fluxes to global matrix
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      for (int d = 0; d < dim; ++d)
        dfdu[q][d] = solution_values[q];
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          viscous_fluxes.add(local_dof_indices[i],
                             local_dof_indices[j],
                             dfdu[q] * fe_values.shape_grad(j, q) *
                               fe_values.shape_value(i, q) * fe_values.JxW(q));
        }
      }
    }
  }
}

/** \brief Computes viscous bilinear forms, to be used in the computation of
 *         maximum-principle preserving first order viscosity.
 *
 *         Each element of the resulting matrix, \f$B_{i,j}\f$ is computed as
 *         follows:
 *         \f[
 *            B_{i,j} = \sum_{K\subset S_{i,j}}b_K(\varphi_i,\varphi_j)
 *         \f]
 */
template <int dim>
void ConservationLaw<dim>::compute_viscous_bilinear_forms()
{
  viscous_bilinear_forms = 0; // zero out matrix
  // local dof indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // query cell volume
    double cell_volume = cell->measure();

    // add local bilinear forms to global matrix
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        double b_cell = 0.0;
        if (j == i)
        {
          b_cell = cell_volume;
        }
        else
        {
          b_cell = -1.0 / (dofs_per_cell - 1.0) * cell_volume;
        }
        viscous_bilinear_forms.add(
          local_dof_indices[i], local_dof_indices[j], b_cell);
      }
    }
  }
}

/**
 * \brief Computes entropy viscosity for each cell.
 *
 * \param[in] dt time step size
 */
template <int dim>
void ConservationLaw<dim>::update_entropy_viscosities(const double & dt)
{
  // compute normalization constant for entropy viscosity
  const double entropy_normalization =
    compute_entropy_normalization(new_solution);

  // get tuning parameters
  const double entropy_residual_coefficient = parameters.entropy_viscosity_coef;
  const double jump_coefficient = parameters.jump_coef;

  // compute entropy viscosity for each cell
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // compute max entropy residual on cell
    const double max_entropy_residual =
      compute_max_entropy_residual(new_solution, old_solution, dt, cell);
    const double max_entropy_jump = compute_max_entropy_jump(new_solution, cell);

    // compute entropy viscosity
    double aux = std::pow(cell_diameter[cell], 2) / entropy_normalization;
    entropy_viscosity[cell] =
      aux * (entropy_residual_coefficient * max_entropy_residual +
             jump_coefficient * max_entropy_jump);
  }
}

/**
 * Computes the domain-averaged entropy and the max entropy deviation in the
 * domain.
 *
 * \param[in] solution solution with which to calculate entropy for the
 *            normalization constant
 *
 * \return normalization constant for entropy viscosity
 */
template <int dim>
double ConservationLaw<dim>::compute_entropy_normalization(
  const Vector<double> & solution) const
{
  FEValues<dim> fe_values(fe, cell_quadrature, update_values | update_JxW_values);

  Vector<double> entropy(n_q_points_cell);

  // compute domain-averaged entropy
  //--------------------------------
  double domain_integral_entropy = 0.0;

  // loop over cells
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // compute entropy
    compute_entropy(solution, fe_values, entropy);

    // loop over quadrature points
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      // add contribution of quadrature point to entropy integral
      domain_integral_entropy += entropy[q] * fe_values.JxW(q);
    }
  }
  // domain-averaged entropy_values
  const double domain_averaged_entropy = domain_integral_entropy / domain_volume;

  // compute max deviation of entropy_values from domain-averaged entropy
  //--------------------------------------------------------------
  double normalization_constant = 0.0;

  // loop over cells
  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // compute entropy
    compute_entropy(solution, fe_values, entropy);

    // loop over quadrature points
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      // update max deviation for the normalization constant
      normalization_constant = std::max(
        normalization_constant, std::abs(entropy[q] - domain_averaged_entropy));
    }
  }

  // guard against division by zero
  if (normalization_constant == 0.0)
    normalization_constant = 1.0;

  return normalization_constant;
}

/**
 * \brief Computes the max entropy residual in a cell.
 *
 * \param[in] new_solution new solution
 * \param[in] old_solution old solution
 * \param[in] dt time step size
 * \param[in] cell cell iterator
 *
 * \return max entropy residual in a cell
 */
template <int dim>
double ConservationLaw<dim>::compute_max_entropy_residual(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const cell_iterator & cell) const
{
  // FE values
  FEValues<dim> fe_values(fe, cell_quadrature, update_values | update_gradients);
  fe_values.reinit(cell);

  Vector<double> entropy_new(n_q_points_cell);
  Vector<double> entropy_old(n_q_points_cell);
  Vector<double> divergence_entropy_flux(n_q_points_cell);
  Vector<double> entropy_residual(n_q_points_cell);

  // compute entropy of current and old solutions
  compute_entropy(new_solution, fe_values, entropy_new);
  compute_entropy(old_solution, fe_values, entropy_old);
  compute_divergence_entropy_flux(
    new_solution, fe_values, divergence_entropy_flux);

  // compute entropy residual at each quadrature point on cell
  double max_entropy_residual = 0.0;
  for (unsigned int q = 0; q < n_q_points_cell; ++q)
  {
    // compute entropy residual
    double dsdt = (entropy_new[q] - entropy_old[q]) / dt;
    entropy_residual[q] = std::abs(dsdt + divergence_entropy_flux[q]);

    // update maximum entropy residual
    max_entropy_residual = std::max(max_entropy_residual, entropy_residual[q]);
  }

  return max_entropy_residual;
}

/**
 * \brief Computes the max entropy jump in a cell.
 *
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 */
template <int dim>
double ConservationLaw<dim>::compute_max_entropy_jump(
  const Vector<double> & solution, const cell_iterator & cell) const
{
  FEFaceValues<dim> fe_values_face(fe,
                                   face_quadrature,
                                   update_values | update_gradients |
                                     update_normal_vectors);

  std::vector<Tensor<1, dim>> gradients_face(n_q_points_face);
  std::vector<Tensor<1, dim>> gradients_face_neighbor(n_q_points_face);
  std::vector<Tensor<1, dim>> normal_vectors(n_q_points_face);
  Vector<double> entropy(n_q_points_face);

  double max_jump_in_cell = 0.0;

  // loop over faces
  for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
  {
    // determine if face is interior
    typename DoFHandler<dim>::face_iterator face = cell->face(iface);
    if (face->at_boundary() == false)
    {
      // get gradients from this cell
      fe_values_face.reinit(cell, iface);
      fe_values_face.get_function_gradients(solution, gradients_face);

      // compute entropy at each quadrature point on face
      compute_entropy(solution, fe_values_face, entropy);

      // get normal vectors
      normal_vectors = fe_values_face.get_all_normal_vectors();

      // get iterator to neighboring cell and face index of this face
      Assert(cell->neighbor(iface).state() == IteratorState::valid,
             ExcInternalError());
      typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(iface);
      const unsigned int ineighbor = cell->neighbor_of_neighbor(iface);
      Assert(ineighbor < faces_per_cell, ExcInternalError());

      // get gradients from neighboring cell
      fe_values_face.reinit(neighbor, ineighbor);
      fe_values_face.get_function_gradients(solution, gradients_face_neighbor);

      // loop over face quadrature points to determine max jump on face
      double max_jump_on_face = 0.0;
      for (unsigned int q = 0; q < n_q_points_face; ++q)
      {
        // compute difference in gradients across face
        gradients_face[q] -= gradients_face_neighbor[q];
        double jump_on_face =
          std::abs((gradients_face[q] * normal_vectors[q]) * entropy[q]);
        max_jump_on_face = std::max(max_jump_on_face, jump_on_face);
      }

      // update max jump in cell
      max_jump_in_cell = std::max(max_jump_in_cell, max_jump_on_face);
    }
  }

  return max_jump_in_cell;
}

/** \brief Computes the negative of the transient residual and stores in
 *         system_rhs
 *  \param i current stage of Runge-Kutta step
 *  \param dt current time step size
 */
/*
template <int dim>
void ConservationLaw<dim>::compute_tr_residual(unsigned int i, double dt)
{
   // the negative of the transient residual is
   // M*(y_current - y_old) - sum{a_ij*dt*f_j}_j=1..i

   // use solution_step to temporarily hold the quantity (y_current - y_old)
   solution_step = new_solution;
   solution_step.add(-1.0,old_solution);

   // store M*(y_current - y_old) in system_rhs
   lumped_mass_matrix.vmult(system_rhs,solution_step);

   // subtract sum{a_ij*dt*f_j}_j=1..i
   for (unsigned int j = 0; j < i; ++j)
      system_rhs.add(-rk.a[i][j]*dt, rk.f[j]);
}
*/

/** \brief Checks that there are no NaNs in the solution vector
 *
 *         The NaN check is performed by comparing a value to itself
 *         since an equality comparison with NaN always returns false.
 */
template <int dim>
void ConservationLaw<dim>::check_nan()
{
  unsigned int n = n_dofs;
  for (unsigned int i = 0; i < n; ++i)
    Assert(new_solution(i) == new_solution(i), ExcNaNEncountered());
}

/** \brief Gets the values and indices of nonzero elements in a sparse matrix.
 *  \param [in] matrix sparse matrix whose row will be retrieved
 *  \param [in] i index of row to be retrieved
 *  \param [out] row_values vector of values of nonzero entries of row i
 *  \param [out] row_indices vector of indices of nonzero entries of row i
 *  \param [out] n_col number of nonzero entries of row i
 */
template <int dim>
void ConservationLaw<dim>::get_matrix_row(const SparseMatrix<double> & matrix,
                                          const unsigned int & i,
                                          std::vector<double> & row_values,
                                          std::vector<unsigned int> & row_indices,
                                          unsigned int & n_col)
{
  // get first and one-past-last iterator for row
  SparseMatrix<double>::const_iterator matrix_iterator = matrix.begin(i);
  SparseMatrix<double>::const_iterator matrix_iterator_end = matrix.end(i);

  // compute number of entries in row and then allocate memory
  n_col = matrix_iterator_end - matrix_iterator;
  row_values.reserve(n_col);
  row_indices.reserve(n_col);

  // loop over columns in row
  for (; matrix_iterator != matrix_iterator_end; ++matrix_iterator)
  {
    row_values.push_back(matrix_iterator->value());
    row_indices.push_back(matrix_iterator->column());
  }
}

/**
 * \brief Checks that the local discrete max principle is satisfied.
 * \param[in] n time step index, used for reporting violations
 */
template <int dim>
bool ConservationLaw<dim>::check_DMP(const unsigned int & n) const
{
  // check that each dof value is bounded by its neighbors
  bool DMP_satisfied = true;
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // check if dof is a Dirichlet node or not
    if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) ==
        dirichlet_nodes.end())
    {
      double value_i = new_solution(i);
      if (value_i < min_values(i))
      {
        DMP_satisfied = false;
        std::cout << "Max principle violated at time step " << n << " with dof "
                  << i << ": " << value_i << " < " << min_values(i) << std::endl;
      }
      if (value_i > max_values(i))
      {
        DMP_satisfied = false;
        std::cout << "Max principle violated at time step " << n << " with dof "
                  << i << ": " << value_i << " > " << max_values(i) << std::endl;
      }
    }
  }

  return DMP_satisfied;
}

/** \brief Computes min and max quantities for max principle
 */
template <int dim>
void ConservationLaw<dim>::compute_max_principle_quantities()
{
  // initialize min and max values
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    max_values(i) = old_solution(i);
    min_values(i) = old_solution(i);
  }

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells
  cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // find min and max values on cell - start by initializing to arbitrary cell
    double max_cell = old_solution(local_dof_indices[0]);
    double min_cell = old_solution(local_dof_indices[0]);
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      double value_j = old_solution(local_dof_indices[j]);
      max_cell = std::max(max_cell, value_j);
      min_cell = std::min(min_cell, value_j);
    }

    // update the max and min values of neighborhood of each dof
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      unsigned int i = local_dof_indices[j]; // global index
      max_values(i) = std::max(max_values(i), max_cell);
      min_values(i) = std::min(min_values(i), min_cell);
    }
  }
}

/**
 * \brief Gets a list of dofs subject to Dirichlet boundary conditions.
 */
template <int dim>
void ConservationLaw<dim>::get_dirichlet_nodes()
{
  // get map of Dirichlet dof indices to Dirichlet values
  std::map<unsigned int, double> boundary_values;

  // clear Dirichlet nodes vector
  dirichlet_nodes.clear();

  if (boundary_conditions_type == "dirichlet")
    // loop over boundary IDs
    for (unsigned int boundary = 0; boundary < n_boundaries; ++boundary)
      // loop over components
      for (unsigned int component = 0; component < n_components; ++component)
      {
        // create mask to prevent function from being applied to all components
        std::vector<bool> component_mask(n_components, false);
        component_mask[component] = true;

        // get boundary values map using interpolate_boundary_values function
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 boundary,
                                                 ZeroFunction<dim>(n_components),
                                                 boundary_values,
                                                 component_mask);

        // extract dof indices from map
        for (std::map<unsigned int, double>::iterator it =
               boundary_values.begin();
             it != boundary_values.end();
             ++it)
          dirichlet_nodes.push_back(it->first);
      }
}

/**
 * \brief Outputs results.
 */
template <int dim>
void ConservationLaw<dim>::output_results(
  PostProcessor<dim> & postprocessor) const
{
  // call post-processor
  postprocessor.output_results(new_solution, dof_handler, triangulation);
}
