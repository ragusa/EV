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
    cout1(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
    timer(mpi_communicator,
          cout,
          TimerOutput::summary,
          TimerOutput::cpu_and_wall_times),
    triangulation(mpi_communicator,
                  typename Triangulation<dim>::MeshSmoothing(
                    Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
#else
    cout1(std::cout, params.verbosity_level >= 1),
    cout2(std::cout, params.verbosity_level >= 2),
    timer(cout2, TimerOutput::summary, TimerOutput::cpu_and_wall_times),
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
    cout1 << std::endl
          << "Cycle " << cycle << " of " << parameters.n_refinement_cycles - 1
          << ":" << std::endl;
    cout1 << "  Number of active cells: " << triangulation.n_active_cells()
          << std::endl;
    cout1 << "  Number of degrees of freedom: " << n_dofs << std::endl;

    // interpolate the initial conditions to the grid, and apply constraints
    VectorTools::interpolate(
      dof_handler, initial_conditions_function, new_solution);
    constraints.distribute(new_solution);

    // set old solutions to the current solution
    old_solution = new_solution;
    older_solution = new_solution;

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

  cout2 << "Initializing system..." << std::endl;

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
    dirichlet_function.resize(n_dirichlet_boundaries);
    for (unsigned int boundary = 0; boundary < n_dirichlet_boundaries; ++boundary)
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

  // determine low-order viscosity and artificial diffusion
  switch (parameters.low_order_scheme)
  {
    case LowOrderScheme::constant:
      low_order_viscosity_type = ViscosityType::constant;
      low_order_diffusion_type = DiffusionType::laplacian;
      break;
    case LowOrderScheme::standard:
      low_order_viscosity_type = ViscosityType::low;
      low_order_diffusion_type = DiffusionType::laplacian;
      break;
    case LowOrderScheme::dmp:
      Assert(false, ExcNotImplemented());
      low_order_viscosity_type = ViscosityType::DMP;
      low_order_diffusion_type = DiffusionType::graphtheoretic;
      break;
    case LowOrderScheme::di_visc:
      low_order_viscosity_type = ViscosityType::DI;
      low_order_diffusion_type = DiffusionType::graphtheoretic;
      break;
    case LowOrderScheme::di_diff:
      low_order_viscosity_type = ViscosityType::none;
      low_order_diffusion_type = DiffusionType::DI;
      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // determine high-order viscosity and artificial diffusion
  switch (parameters.high_order_scheme)
  {
    case HighOrderScheme::galerkin:
      entropy_viscosity_type = ViscosityType::none;
      high_order_viscosity_type = ViscosityType::none;
      high_order_diffusion_type = DiffusionType::none;
      break;
    case HighOrderScheme::entropy_visc:
      entropy_viscosity_type = ViscosityType::entropy;
      high_order_viscosity_type = ViscosityType::high;
      high_order_diffusion_type = DiffusionType::laplacian;
      break;
    case HighOrderScheme::entropy_diff:
      Assert(false, ExcNotImplemented());
      entropy_viscosity_type = ViscosityType::none;
      high_order_viscosity_type = ViscosityType::none;
      high_order_diffusion_type = DiffusionType::entropy;
      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // set unneeded viscosities and diffusion to none
  if (parameters.scheme == Scheme::low)
  {
    entropy_viscosity_type = ViscosityType::none;
    high_order_viscosity_type = ViscosityType::none;
    high_order_diffusion_type = DiffusionType::none;
  }
  if (parameters.scheme == Scheme::high)
  {
    if (parameters.high_order_scheme == HighOrderScheme::galerkin)
    {
      low_order_viscosity_type = ViscosityType::none;
      low_order_diffusion_type = DiffusionType::none;
    }
  }

  // create gradient matrix
  gradient_matrix = std::make_shared<GradientMatrix<dim>>(
    triangulation, cell_quadrature, n_components);

  // create star state
  star_state = create_star_state();
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

/**
 * \brief Adaptively refines mesh based on estimated error per cell.
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

/**
 * \brief Sets up the system before solving.
 *
 *        This function is to be applied after each refinement. It
 *        allocates memory, sets up constraints, makes the sparsity pattern,
 *        and reinitializes the system matrix with the sparsity pattern.
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

  // perform non-standard setup required by derived classes
  perform_nonstandard_setup();

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
    for (unsigned int boundary = 0; boundary < n_dirichlet_boundaries; ++boundary)
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
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
  lumped_mass_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
  low_order_diffusion_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
  high_order_diffusion_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
#else
  // reinitialize matrices with sparsity pattern
  consistent_mass_matrix.reinit(unconstrained_sparsity_pattern);
  lumped_mass_matrix.reinit(unconstrained_sparsity_pattern);
  low_order_diffusion_matrix.reinit(unconstrained_sparsity_pattern);
  high_order_diffusion_matrix.reinit(unconstrained_sparsity_pattern);
#endif

  // assemble mass matrix
  assemble_mass_matrix();

// reinitialize vectors
#ifdef IS_PARALLEL
  old_solution.reinit(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  older_solution.reinit(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  old_stage_solution.reinit(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  new_solution.reinit(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  solution_step.reinit(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  ss_flux.reinit(locally_owned_dofs, mpi_communicator);
  ss_reaction.reinit(locally_owned_dofs, mpi_communicator);
  ss_rhs.reinit(locally_owned_dofs, mpi_communicator);
  estimated_error_per_cell.reinit(triangulation.n_active_cells());
#else
  old_solution.reinit(n_dofs);
  older_solution.reinit(n_dofs);
  old_stage_solution.reinit(n_dofs);
  new_solution.reinit(n_dofs);
  solution_step.reinit(n_dofs);
  system_rhs.reinit(n_dofs);
  ss_flux.reinit(n_dofs);
  ss_reaction.reinit(n_dofs);
  ss_rhs.reinit(n_dofs);
  estimated_error_per_cell.reinit(triangulation.n_active_cells());
#endif

  // reinitialize objects
  star_state->reinitialize();

  // create low-order viscosity
  switch (low_order_viscosity_type)
  {
    // no viscosity
    case ViscosityType::none:
    {
      low_order_viscosity =
        std::make_shared<ConstantViscosity<dim>>(0.0, dof_handler);
      break;
    }
    // constant viscosity
    case ViscosityType::constant:
    {
      low_order_viscosity = std::make_shared<ConstantViscosity<dim>>(
        parameters.constant_viscosity_value, dof_handler);
      break;
    }
    // low-order viscosity
    case ViscosityType::low:
    {
      // ensure diffusion type is compatible
      Assert(low_order_diffusion_type == DiffusionType::none ||
               low_order_diffusion_type == DiffusionType::laplacian,
             ExcInvalidDiffusionType());

      // create viscosity multiplier
      auto viscosity_multiplier = create_viscosity_multiplier();

      // create low-order viscosity
      low_order_viscosity = std::make_shared<LowOrderViscosity<dim>>(
        parameters.first_order_viscosity_coef,
        cell_diameter,
        max_flux_speed_cell,
        fe,
        dof_handler,
        cell_quadrature,
        viscosity_multiplier);

      break;
    }
    // domain-invariant low-order viscosity
    case ViscosityType::DI:
    {
      // ensure diffusion type is compatible
      Assert(low_order_diffusion_type == DiffusionType::none ||
               low_order_diffusion_type == DiffusionType::graphtheoretic,
             ExcInvalidDiffusionType());

      // create viscosity multiplier
      auto viscosity_multiplier = create_viscosity_multiplier();

      // create max wave speed
      auto max_wave_speed = create_max_wave_speed();

      // create domain-invariant viscosity
      low_order_viscosity =
        std::make_shared<DomainInvariantViscosity<dim>>(max_wave_speed,
                                                        gradient_matrix,
                                                        fe,
                                                        dof_handler,
                                                        triangulation,
                                                        cell_quadrature,
                                                        n_components,
                                                        viscosity_multiplier);

      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }

  // create entropy viscosity
  switch (entropy_viscosity_type)
  {
    // no viscosity
    case ViscosityType::none:
    {
      entropy_viscosity =
        std::make_shared<ConstantViscosity<dim>>(0.0, dof_handler);
      break;
    }
    // entropy viscosity
    case ViscosityType::entropy:
    {
      // only implemented for Laplacian diffusion for now
      Assert(high_order_diffusion_type == DiffusionType::laplacian,
             ExcInvalidDiffusionType());

      // create viscosity multiplier
      auto viscosity_multiplier = create_viscosity_multiplier();

      // create entropy
      auto entropy = create_entropy();

      // create entropy viscosity
      entropy_viscosity =
        std::make_shared<EntropyViscosity<dim>>(parameters,
                                                entropy,
                                                cell_diameter,
                                                fe,
                                                dof_handler,
                                                cell_quadrature,
                                                face_quadrature);
      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }

  // create high-order viscosity
  switch (high_order_viscosity_type)
  {
    // no viscosity
    case ViscosityType::none:
    {
      high_order_viscosity =
        std::make_shared<ConstantViscosity<dim>>(0.0, dof_handler);
      break;
    }
    // high-order viscosity
    case ViscosityType::high:
    {
      // ensure that there is a low-order viscosity
      Assert(low_order_viscosity_type != ViscosityType::none, ExcInvalidState());

      // ensure that there is an entropy viscosity
      Assert(entropy_viscosity_type != ViscosityType::none, ExcInvalidState());

      // ensure that low-order viscosity has the same diffusion type
      // since entropy and low-order viscosities are compared to create the
      // high-order viscosity
      Assert(low_order_diffusion_type == high_order_diffusion_type,
             ExcInvalidState());

      // create high-order viscosity
      high_order_viscosity = std::make_shared<HighOrderViscosity<dim>>(
        low_order_viscosity,
        entropy_viscosity,
        parameters.use_low_order_viscosity_for_first_time_step,
        dof_handler);

      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }

  // create low-order and high-order artificial diffusion
  low_order_diffusion = create_artificial_diffusion(low_order_diffusion_type);
  high_order_diffusion = create_artificial_diffusion(high_order_diffusion_type);

  // get list of Dirichlet DoF indices
  get_dirichlet_dof_indices(dirichlet_dof_indices);
}

/**
 * \brief Creates artificial diffusion object
 *
 * \param[in] diffusion_type type of artificial diffusion
 */
template <int dim>
std::shared_ptr<ArtificialDiffusion<dim>> ConservationLaw<
  dim>::create_artificial_diffusion(const DiffusionType & diffusion_type)
{
  std::shared_ptr<ArtificialDiffusion<dim>> artificial_diffusion;

  // create artificial diffusion
  switch (diffusion_type)
  {
    // no diffusion
    case DiffusionType::none:
    {
      artificial_diffusion = std::make_shared<NoDiffusion<dim>>();
      break;
    }
    // Laplacian diffusion
    case DiffusionType::laplacian:
    {
      // get extractors
      std::vector<FEValuesExtractors::Scalar> scalar_extractors;
      std::vector<FEValuesExtractors::Vector> vector_extractors;
      get_fe_extractors(scalar_extractors, vector_extractors);

      artificial_diffusion =
        std::make_shared<LaplacianDiffusion<dim>>(scalar_extractors,
                                                  vector_extractors,
                                                  dof_handler,
                                                  fe,
                                                  cell_quadrature,
                                                  dofs_per_cell);
      break;
    }
    // graph-theoretic diffusion
    case DiffusionType::graphtheoretic:
    {
      artificial_diffusion = std::make_shared<GraphTheoreticDiffusion<dim>>(
        dof_handler, dofs_per_cell, n_components);
      break;
    }
    // domain-invariant diffusion
    case DiffusionType::DI:
    {
      // create max wave speed
      auto max_wave_speed = create_max_wave_speed();

      artificial_diffusion = std::make_shared<DomainInvariantDiffusion<dim>>(
        max_wave_speed, gradient_matrix, n_components, n_dofs);
      break;
    }
    // else
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }

  return artificial_diffusion;
}

/**
 * \brief Updates the cell sizes map and minimum cell size.
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

/**
 * \brief Assembles the mass matrices and applies constraints.
 */
template <int dim>
void ConservationLaw<dim>::assemble_mass_matrix()
{
  // assemble lumped mass matrices (constrained and unconstrained)
  //----------------------------
  assemble_lumped_mass_matrix();

  // assemble consistent mass matrices (constrained and unconstrained)
  //----------------------------
  Function<dim> * dummy_function = 0;
  /*
  MatrixTools::create_mass_matrix(dof_handler,
                                  cell_quadrature,
                                  consistent_mass_matrix_constrained,
                                  dummy_function,
                                  constraints);
  */
  MatrixTools::create_mass_matrix(
    dof_handler, cell_quadrature, consistent_mass_matrix, dummy_function);

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
 *   + \Delta t\sum\limits^{i-1}_{j=1}a_{i,j} \mathbf{r}_j
 * \f]
 *
 * \param[in,out] postprocessor post-processor for outputting transient
 */
template <int dim>
void ConservationLaw<dim>::solve_runge_kutta(PostProcessor<dim> & postprocessor)
{
  // compute steady-state reaction vector, which is not time-dependent
  compute_ss_reaction(ss_reaction);

  // create linear solver
  LinearSolver<dim> linear_solver(
    parameters.linear_solver, constraints, dof_handler, dirichlet_function);

  // create SSPRK time integrator
  SSPRKTimeIntegrator<dim> ssprk(parameters.time_discretization,
                                 n_dofs,
                                 linear_solver,
                                 unconstrained_sparsity_pattern);

  // create FCT if applicable
  std::shared_ptr<FCT<dim>> fct;
  if (parameters.scheme == Scheme::fct)
    fct = std::make_shared<FCT<dim>>(parameters,
                                     dof_handler,
                                     triangulation,
                                     lumped_mass_matrix,
                                     consistent_mass_matrix,
                                     star_state,
                                     linear_solver,
                                     unconstrained_sparsity_pattern,
                                     dirichlet_dof_indices,
                                     n_components,
                                     dofs_per_cell,
                                     component_names);

  // initialize old time, time index, and transient flag
  old_time = 0.0;
  double old_dt = 0.0;
  unsigned int n = 1;
  bool in_transient = true;

  // start transient
  while (in_transient)
  {
    // update max speed for use in CFL computation
    update_flux_speeds();

    // compute dt
    double dt;
    switch (parameters.time_step_size_method)
    {
      case TimeStepSizeMethod::constant:
        dt = parameters.time_step_size;
        break;
      case TimeStepSizeMethod::cfl:
        dt = compute_dt_from_cfl_condition();
        break;
      default:
        Assert(false, ExcNotImplemented());
        break;
    }
    // check end of transient and shorten last time step if necessary
    if ((old_time + dt) >= end_time)
    {
      dt = end_time - old_time;
      in_transient = false;
    }

    // log the time step size for computing average time step size for
    // convergence table
    postprocessor.log_time_step_size(dt);

    // compute new time
    new_time = old_time + dt;

    // compute CFL number for printing
    const double cfl = compute_cfl_number(dt);
    const bool cfl_is_violated = cfl > parameters.cfl + 1.0e-15;

    // print CFL in red if it violates CFL condition
    if (cfl_is_violated)
      cout1 << std::fixed << std::setprecision(2) << "  time step " << n
            << ": t = \x1b[1;34m" << new_time << "\x1b[0m, CFL = \x1b[1;31m"
            << cfl << "\x1b[0m" << std::endl;
    else
      cout1 << std::fixed << std::setprecision(2) << "  time step " << n
            << ": t = \x1b[1;34m" << new_time << "\x1b[0m, CFL = " << cfl
            << std::endl;

    // initialize SSPRK time step
    ssprk.initialize_time_step(old_solution, dt);

    // loop over SSPRK stages to compute new solution
    for (unsigned int i = 0; i < ssprk.n_stages; ++i)
    {
      // print stage index
      cout1 << "    stage " << i + 1 << " of " << ssprk.n_stages << std::endl;

      // get stage time
      double t_stage = ssprk.get_stage_time();

      // compute steady-state rhs vector
      compute_ss_rhs(t_stage, ss_rhs);

      // compute steady-state flux vector
      compute_ss_flux(dt, old_solution, ss_flux);

      // determine old stage solution and dt (needed for entropy viscosity)
      double old_stage_dt;
      if (i == 0)
      {
        old_stage_dt = old_dt;
      }
      else
      {
        ssprk.get_stage_solution(i - 1, older_solution);
        old_stage_dt = dt;
      }

      // get old stage solution
      ssprk.get_stage_solution(i, old_stage_solution);

      // advance by an SSPRK step
      switch (parameters.scheme)
      {
        case Scheme::low:
          // update star states
          // star_state->compute_star_states(old_stage_solution);

          low_order_viscosity->update(
            old_stage_solution, older_solution, old_stage_dt, n);
          low_order_diffusion->compute_diffusion_matrix(
            old_stage_solution, low_order_viscosity, low_order_diffusion_matrix);
          ssprk.step(lumped_mass_matrix,
                     ss_flux,
                     ss_rhs,
                     low_order_diffusion_matrix,
                     true);
          break;
        case Scheme::high:
          high_order_viscosity->update(
            old_stage_solution, older_solution, old_stage_dt, n);
          high_order_diffusion->compute_diffusion_matrix(
            old_stage_solution,
            high_order_viscosity,
            high_order_diffusion_matrix);
          ssprk.step(consistent_mass_matrix,
                     ss_flux,
                     ss_rhs,
                     high_order_diffusion_matrix,
                     true);
          break;
        case Scheme::fct:
          perform_fct_ssprk_step(dt, old_stage_dt, n, fct, ssprk);
          break;
        default:
          Assert(false, ExcNotImplemented());
          break;
      }

      // save viscosity if it is the first stage
      if (i == 0)
        if (parameters.output_viscosity_transient)
          output_viscosity(postprocessor, true, old_time);
    }

    // retrieve the final solution
    ssprk.get_new_solution(new_solution);

    {
      // start timer for output
      TimerOutput::Scope timer_section(timer, "Output");

      // output solution transient if specified
      postprocessor.output_solution_transient(
        new_solution, new_time, dof_handler, "solution", false);

      // output FCT bounds if specified
      if (fct && parameters.output_fct_bounds)
        fct->output_bounds_transient(postprocessor, new_time);
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
    cout1 << std::scientific << std::setprecision(4)
          << "    Steady-state error = "
          << "\x1b[1;35m" << steady_state_error << "\x1b[0m" << std::endl;
    if (steady_state_error < parameters.steady_state_tolerance)
    {
      in_transient = false;
      cout1 << std::scientific << std::setprecision(4) << "\x1b[1;32m"
            << "  Converged to steady-state."
            << "\x1b[0m" << std::endl;
    }

    // compute error for adaptive mesh refinement
    compute_error_for_refinement();

    // increment transient counter for post-processor
    postprocessor.increment_transient_counter();

    // store old solutions
    older_solution = old_solution;
    old_solution = new_solution;

    // reset old time and increment time step index
    old_time = new_time;
    old_dt = dt;
    n++;
  } // end of time loop

  // output limiter matrix if specified
  if (parameters.scheme == Scheme::fct)
    if (parameters.output_limiter_matrix)
      fct->output_limiter_matrix();
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
  if (low_order_viscosity_type != ViscosityType::none)
  {
    viscosities.push_back(low_order_viscosity);
    viscosity_names.push_back("low_order_viscosity");
  }
  if (entropy_viscosity_type != ViscosityType::none)
  {
    viscosities.push_back(low_order_viscosity);
    viscosity_names.push_back("entropy_viscosity");
  }
  if (high_order_viscosity_type != ViscosityType::none)
  {
    viscosities.push_back(high_order_viscosity);
    viscosity_names.push_back("high_order_viscosity");
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

/**
 * \brief Creates a default viscosity multiplier object, which uses the number
 *        1 as the multiplier.
 *
 * \return pointer to created viscosity multiplier
 */
template <int dim>
std::shared_ptr<ViscosityMultiplier<dim>> ConservationLaw<
  dim>::create_viscosity_multiplier() const
{
  auto viscosity_multiplier = std::make_shared<ViscosityMultiplier<dim>>();
  return viscosity_multiplier;
}

/**
 * \brief Performs an SSPRK step using FCT.
 *
 * \param[in] dt current time step size
 * \param[in] old_stage_dt time step size of previous SSPRK stage
 * \param[in] n time index
 * \param[in] fct FCT
 * \param[in,out] ssprk SSPRK time integrator
 */
template <int dim>
void ConservationLaw<dim>::perform_fct_ssprk_step(
  const double & dt,
  const double & old_stage_dt,
  const unsigned int & n,
  const std::shared_ptr<FCT<dim>> & fct,
  SSPRKTimeIntegrator<dim> & ssprk)
{
  // update star states
  star_state->compute_star_states(old_stage_solution);

  // update low-order diffusion
  low_order_viscosity->update(
    old_stage_solution, older_solution, old_stage_dt, n);
  low_order_diffusion->compute_diffusion_matrix(
    old_stage_solution, low_order_viscosity, low_order_diffusion_matrix);

  // update high-order diffusion
  high_order_viscosity->update(
    old_stage_solution, older_solution, old_stage_dt, n);
  high_order_diffusion->compute_diffusion_matrix(
    old_stage_solution, high_order_viscosity, high_order_diffusion_matrix);

  // compute high-order solution
  ssprk.step(
    consistent_mass_matrix, ss_flux, ss_rhs, high_order_diffusion_matrix, false);

  // get high-order solution
  ssprk.get_intermediate_solution(new_solution);

  // perform FCT
  fct->solve_fct_system(new_solution,
                        old_stage_solution,
                        ss_flux,
                        ss_reaction,
                        ss_rhs,
                        dt,
                        low_order_diffusion_matrix,
                        high_order_diffusion_matrix);

  // set stage solution to be FCT solution for this stage
  ssprk.set_intermediate_solution(new_solution);

  // finish computing stage solution
  ssprk.complete_stage_solution();
}

/**
 * \brief Gets a vector of dof indices subject to Dirichlet boundary conditions
 *
 * \param[out] dirichlet_dof_indices vector of dof indices subject to Dirichlet
 *             boundary conditions
 */
template <int dim>
void ConservationLaw<dim>::get_dirichlet_dof_indices(
  std::vector<unsigned int> & dirichlet_dof_indices)
{
  if (boundary_conditions_type == "dirichlet")
  {
    // get map of Dirichlet dof indices to Dirichlet values
    std::map<unsigned int, double> boundary_values;

    // loop over boundary IDs
    for (unsigned int boundary = 0; boundary < n_dirichlet_boundaries; ++boundary)
      // loop over components
      for (unsigned int component = 0; component < n_components; ++component)
      {
        // mask other components
        std::vector<bool> component_mask(n_components, false);
        component_mask[component] = true;

        // fill boundary_values with boundary values
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 boundary,
                                                 *(dirichlet_function[boundary]),
                                                 boundary_values,
                                                 component_mask);
      }

    // extract dof indices from map
    dirichlet_dof_indices.clear();
    for (std::map<unsigned int, double>::iterator it = boundary_values.begin();
         it != boundary_values.end();
         ++it)
      dirichlet_dof_indices.push_back(it->first);
  }
}
