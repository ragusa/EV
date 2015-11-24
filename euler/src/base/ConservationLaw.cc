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
  :
#ifdef IS_PARALLEL
    mpi_communicator(MPI_COMM_WORLD),
    cout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
    timer(mpi_communicator,
          cout,
          TimerOutput::summary,
          TimerOutput::cpu_and_wall_times),
    triangulation(mpi_communicator,
                  typename Triangulation<dim>::MeshSmoothing(
                    Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
#else
    cout(std::cout, params.verbosity_level > 0),
    timer(cout, TimerOutput::summary, TimerOutput::cpu_and_wall_times),
    triangulation(),
#endif
    parameters(params),
    n_components(params.n_components),
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

  // create data post-processor (for outputting auxiliary variables)
  std::shared_ptr<DataPostprocessor<dim>> aux_postprocessor =
    create_auxiliary_postprocessor();

  // create post-processor object
  PostProcessor<dim> postprocessor(parameters,
                                   end_time,
                                   has_exact_solution,
                                   exact_solution_function,
                                   problem_name,
                                   component_names,
                                   component_interpretations,
                                   triangulation,
                                   aux_postprocessor);

  // loop over adaptive refinement cycles
  for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; ++cycle)
  {
    // set cycle for post-processor
    postprocessor.set_cycle(cycle);

    // refine mesh if not the first cycle
    if (cycle > 0)
      refine_mesh();

    // setup system; to be applied after each refinement
    setup_system();

    // print information
    cout << std::endl
         << "Cycle " << cycle << " of " << parameters.n_refinement_cycles - 1
         << ":" << std::endl;
    cout << "  Number of active cells: " << triangulation.n_active_cells()
         << std::endl;

    // interpolate the initial conditions to the grid, and apply constraints
    VectorTools::interpolate(
      dof_handler, initial_conditions_function, new_solution);
    apply_Dirichlet_BC(0.0);
    constraints.distribute(new_solution);

    // set old solution to the current solution
    old_solution = new_solution;

    // output block
    {
      // start timer for output
      TimerOutput::Scope timer_section(timer, "Output");

      // output initial solution
      postprocessor.output_solution(
        new_solution, 0.0, dof_handler, "solution_initial", false);
      postprocessor.output_solution_transient(
        new_solution, 0.0, dof_handler, "solution", false);
    }

    // solve transient with selected time integrator
    switch (parameters.temporal_integrator)
    {
      case TemporalIntegrator::runge_kutta: // Runge-Kutta
        solve_runge_kutta(postprocessor);
        break;
      default:
        Assert(false, ExcNotImplemented());
        break;
    }

    // evaluate error block
    {
      // start timer for error evaluation
      TimerOutput::Scope timer_section(timer, "Evaluate error");

      // evaluate errors for convergence study
      postprocessor.evaluate_error(new_solution, dof_handler, triangulation);
    }

  } // end of refinement loop

  // output block
  {
    // start timer for output
    TimerOutput::Scope timer_section(timer, "Output");

    // call post-processor to output solution and convergence results
    postprocessor.output_results(new_solution, dof_handler, triangulation);
  }

  // output viscosity if requested
  if (parameters.output_viscosity)
    output_viscosity(postprocessor);

  // print final solution if requested
  if (parameters.print_final_solution)
  {
    std::cout.precision(10);
    std::cout.setf(std::ios::scientific);

    for (unsigned int i = 0; i < n_dofs; ++i)
      std::cout << new_solution[i] << std::endl;
  }
}

/**
 * \brief Initializes system.
 */
template <int dim>
void ConservationLaw<dim>::initialize_system()
{
  // start timer for initialize section
  TimerOutput::Scope timer_section(timer, "Initialize");

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
  if (parameters.temporal_integrator == TemporalIntegrator::runge_kutta)
    initialize_runge_kutta();

  // determine mass matrix to be used based on method
  switch (parameters.viscosity_type)
  {
    case ViscosityType::none:
      mass_matrix = &consistent_mass_matrix;
      break;
    case ViscosityType::constant:
      mass_matrix = &lumped_mass_matrix;
      break;
    case ViscosityType::low:
      mass_matrix = &lumped_mass_matrix;
      break;
    case ViscosityType::DMP_low:
      mass_matrix = &lumped_mass_matrix;
      break;
    case ViscosityType::entropy:
      mass_matrix = &consistent_mass_matrix;
      break;
    default:
      ExcNotImplemented();
      break;
  }
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
    case TemporalDiscretization::ERK1:
      rk.s = 1;
      break;
    case TemporalDiscretization::ERK2:
      rk.s = 2;
      break;
    case TemporalDiscretization::ERK3:
      rk.s = 3;
      break;
    case TemporalDiscretization::ERK4:
      rk.s = 4;
      break;
    case TemporalDiscretization::SDIRK22:
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
    case TemporalDiscretization::ERK1:
      rk.b[0] = 1;
      rk.c[0] = 0;
      rk.is_explicit = true;
      break;
    case TemporalDiscretization::ERK2:
      rk.a[1][0] = 0.5;
      rk.b[0] = 0;
      rk.b[1] = 1;
      rk.c[0] = 0;
      rk.c[1] = 0.5;
      rk.is_explicit = true;
      break;
    case TemporalDiscretization::ERK3:
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
    case TemporalDiscretization::ERK4:
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
    case TemporalDiscretization::SDIRK22:
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

/**
 * \brief Computes error for adaptive mesh refinement for a time
 *        step and adds it to an error sum for all time steps.
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
  // start timer for refinement section
  TimerOutput::Scope timer_section(timer, "Refinement");

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
  // start timer for setup section
  TimerOutput::Scope timer_section(timer, "Setup");

  // clear maps
  cell_diameter.clear();
  max_flux_speed_cell.clear();

  // clear and distribute dofs
  dof_handler.clear();
  dof_handler.distribute_dofs(fe);
  n_dofs = dof_handler.n_dofs();
  cout << "Number of degrees of freedom: " << n_dofs << std::endl;

#ifdef IS_PARALLEL
  // get index set of locally owned DoFs
  locally_owned_dofs = dof_handler.locally_owned_dofs();

  // get index set of locally relevant DoFs
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
#endif

  // update cell sizes and minimum cell size
  update_cell_sizes();

  // make constraints
  constraints.clear();
#ifdef IS_PARALLEL
  // tell constraint matrix the list of locally relevant DoFs
  constraints.reinit(locally_relevant_dofs);
#endif
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

// make dynamic sparsity pattern
#ifdef IS_PARALLEL
  DynamicSparsityPattern dsp_constrained(locally_relevant_dofs);
  DynamicSparsityPattern dsp_unconstrained(locally_relevant_dofs);
#else
  DynamicSparsityPattern dsp_constrained(n_dofs);
  DynamicSparsityPattern dsp_unconstrained(n_dofs);
#endif

  // make sparsity pattern
  DoFTools::make_sparsity_pattern(
    dof_handler, dsp_constrained, constraints, false);
  DoFTools::make_sparsity_pattern(dof_handler, dsp_unconstrained);
  constrained_sparsity_pattern.copy_from(dsp_constrained);
  unconstrained_sparsity_pattern.copy_from(dsp_unconstrained);

#ifdef IS_PARALLEL
  // distribute sparsity pattern
  SparsityTools::distribute_sparsity_pattern(
    dsp_constrained,
    dof_handler.n_locally_owned_dofs_per_processor(),
    mpi_communicator,
    locally_relevant_dofs);
  SparsityTools::distribute_sparsity_pattern(
    dsp_unconstrained,
    dof_handler.n_locally_owned_dofs_per_processor(),
    mpi_communicator,
    locally_relevant_dofs);

  // reinitialize matrices with sparsity pattern
  consistent_mass_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_constrained, mpi_communicator);
  lumped_mass_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_constrained, mpi_communicator);
  system_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_constrained, mpi_communicator);
#else
  // reinitialize matrices with sparsity pattern
  consistent_mass_matrix.reinit(constrained_sparsity_pattern);
  lumped_mass_matrix.reinit(constrained_sparsity_pattern);
  system_matrix.reinit(constrained_sparsity_pattern);
#endif

  // assemble mass matrix
  assemble_mass_matrix();

// reinitialize vectors
#ifdef IS_PARALLEL
  old_solution.reinit(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  new_solution.reinit(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  solution_step.reinit(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  estimated_error_per_cell.reinit(triangulation.n_active_cells());
#else
  old_solution.reinit(n_dofs);
  new_solution.reinit(n_dofs);
  solution_step.reinit(n_dofs);
  system_rhs.reinit(n_dofs);
  estimated_error_per_cell.reinit(triangulation.n_active_cells());
#endif

  // allocate memory for steady-state residual evaluations if Runge-Kutta is
  // used
  if (parameters.temporal_integrator == TemporalIntegrator::runge_kutta)
    for (int i = 0; i < rk.s; ++i)
#ifdef IS_PARALLEL
      rk.f[i].reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
#else
      rk.f[i].reinit(n_dofs);
#endif

  // create viscosity
  switch (parameters.viscosity_type)
  {
    // no viscosity
    case ViscosityType::none:
    {
      auto zero_viscosity_tmp =
        std::make_shared<ConstantViscosity<dim>>(0.0, dof_handler);
      viscosity = zero_viscosity_tmp;
      break;
    }
    // constant viscosity
    case ViscosityType::constant:
    {
      auto constant_viscosity_tmp = std::make_shared<ConstantViscosity<dim>>(
        parameters.constant_viscosity_value, dof_handler);
      viscosity = constant_viscosity_tmp;
      break;
    }
    // low-order viscosity
    case ViscosityType::low:
    {
      auto low_order_viscosity_tmp = std::make_shared<LowOrderViscosity<dim>>(
        parameters.first_order_viscosity_coef,
        cell_diameter,
        max_flux_speed_cell,
        dof_handler);
      low_order_viscosity = low_order_viscosity_tmp;

      viscosity = low_order_viscosity_tmp;

      break;
    }
    // entropy viscosity
    case ViscosityType::entropy:
    {
      // create entropy
      auto entropy = create_entropy();

      // create low-order viscosity
      auto low_order_viscosity_tmp = std::make_shared<LowOrderViscosity<dim>>(
        parameters.first_order_viscosity_coef,
        cell_diameter,
        max_flux_speed_cell,
        dof_handler);
      low_order_viscosity = low_order_viscosity_tmp;

      // create entropy viscosity
      auto entropy_viscosity_tmp =
        std::make_shared<EntropyViscosity<dim>>(parameters,
                                                entropy,
                                                cell_diameter,
                                                fe,
                                                dof_handler,
                                                cell_quadrature,
                                                face_quadrature);
      entropy_viscosity = entropy_viscosity_tmp;

      // create high-order viscosity
      auto high_order_viscosity_tmp = std::make_shared<HighOrderViscosity<dim>>(
        low_order_viscosity,
        entropy_viscosity,
        parameters.use_low_order_viscosity_for_first_time_step,
        dof_handler);

      viscosity = high_order_viscosity_tmp;

      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }

  // create artificial diffusion
  switch (parameters.diffusion_type)
  {
    // no diffusion
    case DiffusionType::none:
    {
      auto tmp_artificial_diffusion = std::make_unique<NoDiffusion<dim>>();
      artificial_diffusion = std::move(tmp_artificial_diffusion);
      break;
    }
    // algebraic diffusion
    case DiffusionType::algebraic:
    {
      auto tmp_artificial_diffusion = std::make_unique<AlgebraicDiffusion<dim>>();
      artificial_diffusion = std::move(tmp_artificial_diffusion);
      break;
    }
    // Laplacian diffusion
    case DiffusionType::laplacian:
    {
      // get extractors
      std::vector<FEValuesExtractors::Scalar> scalar_extractors;
      std::vector<FEValuesExtractors::Vector> vector_extractors;
      get_fe_extractors(scalar_extractors, vector_extractors);

      auto tmp_artificial_diffusion = std::make_unique<LaplacianDiffusion<dim>>(
        scalar_extractors, vector_extractors, n_q_points_cell, dofs_per_cell);
      artificial_diffusion = std::move(tmp_artificial_diffusion);
      break;
    }
    // graph-theoretic diffusion
    case DiffusionType::graphtheoretic:
    {
      ExcNotImplemented();
      break;
    }
    // else
    default:
    {
      ExcNotImplemented();
      break;
    }
  }

  // perform additional setup required by derived classes
  perform_additional_setup();
}

/** \brief Updates the cell sizes map and minimum cell size.
 */
template <int dim>
void ConservationLaw<dim>::update_cell_sizes()
{
  // fill cell size map and find minimum cell size
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
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
 * \param[inout] postprocessor post-processor for outputting transient
 */
template <int dim>
void ConservationLaw<dim>::solve_runge_kutta(PostProcessor<dim> & postprocessor)
{
  old_time = 0.0;
  unsigned int n = 1; // time step index
  double t_end = end_time;
  bool in_transient = true;
  // bool DMP_satisfied = true;
  while (in_transient)
  {
    // update max speed for use in CFL computation
    update_flux_speeds();

    // compute dt
    double dt;
    switch (parameters.time_step_size_method)
    {
      case TimeStepSizeMethod::constant_dt:
        dt = parameters.time_step_size;
        break;
      case TimeStepSizeMethod::cfl_condition:
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

    // compute new time
    new_time = old_time + dt;

    // compute CFL number for printing
    const double cfl = compute_cfl_number(dt);
    const bool cfl_is_violated = cfl > parameters.cfl + 1.0e-15;

    // print CFL in red if it violates CFL condition
    if (cfl_is_violated)
      cout << std::fixed << std::setprecision(2) << "  time step " << n
           << ": t = \x1b[1;34m" << new_time << "\x1b[0m, CFL = \x1b[1;31m" << cfl
           << "\x1b[0m" << std::endl;
    else
      cout << std::fixed << std::setprecision(2) << "  time step " << n
           << ": t = \x1b[1;34m" << new_time << "\x1b[0m, CFL = " << cfl
           << std::endl;

    // compute viscosities
    {
      // start timer for compute viscosities section
      TimerOutput::Scope timer_section(timer, "Compute viscosity");

      // update viscosity
      viscosity->update(new_solution, old_solution, dt, n);
    }

    // output viscosity transient if specified
    if (parameters.output_viscosity_transient)
      output_viscosity(postprocessor, true /* is_transient */, old_time);

    // update old_solution to new_solution for next time step;
    // this is not done at the end of the previous time step because
    // of the time derivative term in update_viscosities() above
    old_solution = new_solution;

    // solve here
    // compute each f_i
    for (int i = 0; i < rk.s; ++i)
    {
      cout << "    stage " << i + 1 << " of " << rk.s << std::endl;

      // compute stage time
      const double stage_time = old_time + rk.c[i] * dt;

      if (rk.is_explicit) // explicit Runge-Kutta
      {
        /* compute RHS of
         *
         * M*Y_i = M*y^n + sum_{j=1}^{i-1} dt*a_{i,j}*f_j
         *
         */
        system_rhs = 0.0;
        mass_matrix->vmult(system_rhs, old_solution);
        for (int j = 0; j < i; ++j)
          system_rhs.add(dt * rk.a[i][j], rk.f[j]);

        // solve system M*Y_i = RHS
        linear_solve(*mass_matrix, system_rhs, new_solution);

        // ordinarily, Dirichlet BC need not be reapplied, but in general,
        // the Dirichlet BC can be time-dependent
        apply_Dirichlet_BC(stage_time);

        // distribute hanging node constraints
        constraints.distribute(new_solution);
      }
      else
      { // implicit Runge-Kutta
        Assert(false, ExcNotImplemented());
      }

      // compute steady-state residual
      {
        // start timer for compute steady-state residual function
        TimerOutput::Scope timer_section(timer, "Compute steady-state residual");

        compute_ss_residual(dt, rk.f[i]);
      }
    }

    /* compute the solution using the computed f_i's:
     *
     * M*y^{n+1} = M*y^n + dt*sum_{i=1}^s b_i*f_i
     */
    // solve M*y^{n+1} = M*y^n + dt*sum_{i=1}^s b_i*f_i
    system_rhs = 0.0;
    // lumped_mass_matrix.vmult(system_rhs, old_solution);
    mass_matrix->vmult(system_rhs, old_solution);
    for (int i = 0; i < rk.s; ++i)
      system_rhs.add(dt * rk.b[i], rk.f[i]);
    // linear_solve(lumped_mass_matrix, system_rhs, new_solution);
    linear_solve(*mass_matrix, system_rhs, new_solution);

    // apply Dirichlet constraints
    apply_Dirichlet_BC(new_time);

    // distribute hanging node constraints
    constraints.distribute(new_solution);

    {
      // start timer for output
      TimerOutput::Scope timer_section(timer, "Output");

      // output solution transient if specified
      postprocessor.output_solution_transient(
        new_solution, new_time, dof_handler, "solution", false);
    }

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
    cout << std::scientific << std::setprecision(4) << "    Steady-state error = "
         << "\x1b[1;35m" << steady_state_error << "\x1b[0m" << std::endl;
    if (steady_state_error < parameters.steady_state_tolerance)
    {
      in_transient = false;
      cout << std::scientific << std::setprecision(4) << "\x1b[1;32m"
           << "  Converged to steady-state."
           << "\x1b[0m" << std::endl;
    }

    /*
        // check that DMP is satisfied at all time steps
        if (parameters.viscosity_type ==
            max_principle)
        {
          bool DMP_satisfied_this_time_step = check_DMP(n);
          DMP_satisfied = DMP_satisfied and DMP_satisfied_this_time_step;
        }
    */

    // compute error for adaptive mesh refinement
    compute_error_for_refinement();

    // reset old time and increment time step index
    old_time = new_time;
    n++;
  } // end of time loop

  /*
    // report if DMP was satisfied at all time steps
    if (parameters.viscosity_type ==
    ViscosityType::DMP_low)
    {
      if (DMP_satisfied)
        cout << "DMP was satisfied at all time steps" << std::endl;
      else
        cout << "DMP was NOT satisfied at all time steps" << std::endl;
    }
  */
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
  // start timer for linear solve section
  TimerOutput::Scope timer_section(timer, "Linear solve");

  switch (parameters.linear_solver)
  {
    // UMFPACK
    case LinearSolverType::direct:
    {
      SparseDirectUMFPACK A_umfpack;
      A_umfpack.initialize(A);
      A_umfpack.vmult(x, b);
      break;
    }
    // GMRes
    case LinearSolverType::gmres:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
    // CG with SSOR
    case LinearSolverType::cg:
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
 * \brief Returns a pointer to an auxiliary post-processor object if there is
 *        one; otherwise returns a nullptr
 *
 * Derived classes define this function if any additional quantities are
 * desired to be output.
 *
 * \return pointer to the auxiliary post-processor object or nullptr
 */
template <int dim>
std::shared_ptr<DataPostprocessor<dim>> ConservationLaw<
  dim>::create_auxiliary_postprocessor() const
{
  return nullptr;
}

/**
 * \brief Outputs viscosities.
 *
 * \param[in] postprocessor postprocessor
 * \param[in] is_transient flag signalling that this is to output as an
 *            item in a transient
 * \param[in] time current time value, which is used if this is output as an
 *            item in a transient
 */
template <int dim>
void ConservationLaw<dim>::output_viscosity(PostProcessor<dim> & postprocessor,
                                            const bool & is_transient,
                                            const double & time)
{
  // start timer for output
  TimerOutput::Scope timer_section(timer, "Output");

  // create vector of pointers to cell maps of viscosity
  std::vector<std::shared_ptr<Viscosity<dim>>> viscosities;
  std::vector<std::string> viscosity_names;

  // output final viscosities if non-constant viscosity used
  switch (parameters.viscosity_type)
  {
    case ViscosityType::none:
      break;
    case ViscosityType::constant:
      break;
    case ViscosityType::low:
      viscosities.push_back(low_order_viscosity);
      viscosity_names.push_back("low_order_viscosity");
      break;
    case ViscosityType::entropy:
      // low-order viscosity
      viscosities.push_back(low_order_viscosity);
      viscosity_names.push_back("low_order_viscosity");
      // entropy viscosity
      viscosities.push_back(entropy_viscosity);
      viscosity_names.push_back("entropy_viscosity");
      // high-order viscosity
      viscosities.push_back(viscosity);
      viscosity_names.push_back("high_order_viscosity");

      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // output the  viscosities
  if (viscosities.size() > 0)
    if (is_transient)
      postprocessor.output_viscosity_transient(
        viscosities, viscosity_names, time, dof_handler);
    else
      postprocessor.output_viscosity(
        viscosities, viscosity_names, time, dof_handler);
}

/**
 * \brief Checks that there are no NaNs in the solution vector
 *
 * The NaN check is performed by comparing a value to itself
 * since an equality comparison with NaN always returns false.
 */
template <int dim>
void ConservationLaw<dim>::check_nan()
{
  unsigned int n = n_dofs;
  for (unsigned int i = 0; i < n; ++i)
    Assert(new_solution(i) == new_solution(i), ExcNaNEncountered());
}
