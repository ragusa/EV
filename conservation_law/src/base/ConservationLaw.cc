/**
 * \file ConservationLaw.cc
 * \brief Provides the function definitions for the ConservationLaw class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] params conservation law parameters
 * \param[in] n_components_  number of components for conservation law
 * \param[in] is_linear_  flag that conservation law is linear
 */
template <int dim>
ConservationLaw<dim>::ConservationLaw(const RunParameters & params,
                                      const unsigned int & n_components_,
                                      const bool & is_linear_)
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
    n_components(n_components_),
    are_star_states(!is_linear_),
    fe(FE_Q<dim>(params.degree), n_components_),
    dofs_per_cell(fe.dofs_per_cell),
    faces_per_cell(GeometryInfo<dim>::faces_per_cell),
    dof_handler(triangulation),
    n_q_points_per_dim(params.n_quadrature_points),
    cell_quadrature(n_q_points_per_dim),
    face_quadrature(n_q_points_per_dim),
    n_q_points_cell(cell_quadrature.size()),
    n_q_points_face(face_quadrature.size()),
    is_scalar(n_components_ == 1),
    is_linear(is_linear_),
    need_to_compute_inviscid_ss_matrix(false)
{
  // assert that an explicit time integrator is specified
  AssertThrow(params.temporal_discretization ==
                TemporalDiscretizationClassification::ssprk,
              ExcNotImplemented());
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

  // compute end time
  double end_time;
  if (parameters.use_default_end_time &&
      problem_base_parameters->has_default_end_time)
    end_time = problem_base_parameters->default_end_time;
  else
    end_time = parameters.end_time;

  // create post-processor object
  PostProcessor<dim> postprocessor(
    parameters,
    n_components,
    end_time,
    problem_base_parameters->has_exact_solution,
    problem_base_parameters->exact_solution_function,
    problem_base_parameters->problem_name,
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
    problem_base_parameters->initial_conditions_function.set_time(0.0);
    VectorTools::interpolate(dof_handler,
                             problem_base_parameters->initial_conditions_function,
                             new_solution);
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

    // solve transient with selected time discretization
    switch (parameters.temporal_discretization)
    {
      case TemporalDiscretizationClassification::ssprk: // SSPRK
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

  // perform initial mesh refinements
  triangulation.refine_global(parameters.initial_refinement_level);

  // determine low-order viscosity and artificial diffusion
  switch (parameters.low_order_scheme)
  {
    case LowOrderScheme::constant:
      low_order_viscosity_type = ViscosityType::constant;
      low_order_diffusion_type = DiffusionType::laplacian;
      break;
    case LowOrderScheme::lax:
      low_order_viscosity_type = ViscosityType::lax;
      low_order_diffusion_type = DiffusionType::laplacian;
      break;
    case LowOrderScheme::dmp:
      // assert that system is scalar
      AssertThrow(is_scalar, ExcInvalidState());

      // DMP low-order viscosity requires computation of inviscid ss matrix
      need_to_compute_inviscid_ss_matrix = true;

      // issue warning if an inappropriate time step size was chosen
      if (parameters.time_step_size_option != TimeStepSizeOption::constant &&
          parameters.time_step_size_option != TimeStepSizeOption::cfl_dmp)
        utilities::issue_warning(
          "The DMP low-order method was chosen for the low-order scheme, \n"
          "but the DMP CFL time step size method was not used. \n"
          "The time step sizes used may violate the CFL condition.");

      low_order_viscosity_type = ViscosityType::DMP;
      low_order_diffusion_type = DiffusionType::graphtheoretic;
      break;
    case LowOrderScheme::di_visc:
      // issue warning if an inappropriate time step size was chosen
      if (parameters.time_step_size_option != TimeStepSizeOption::constant &&
          parameters.time_step_size_option != TimeStepSizeOption::cfl_di)
        utilities::issue_warning(
          "The DI low-order method was chosen for the low-order scheme, \n"
          "but the DI CFL time step size method was not used. \n"
          "The time step sizes used may violate the CFL condition");

      low_order_viscosity_type = ViscosityType::DI;
      low_order_diffusion_type = DiffusionType::graphtheoretic;
      break;
    case LowOrderScheme::di_diff:
      // issue warning if an inappropriate time step size was chosen
      if (parameters.time_step_size_option != TimeStepSizeOption::constant &&
          parameters.time_step_size_option != TimeStepSizeOption::cfl_di)
        utilities::issue_warning(
          "The DI low-order method was chosen for the low-order scheme, \n"
          "but the DI CFL time step size method was not used. \n"
          "The time step sizes used may violate the CFL condition");

      low_order_viscosity_type = ViscosityType::none;
      low_order_diffusion_type = DiffusionType::DI;
      break;
    default:
      AssertThrow(false, ExcNotImplemented());
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
      high_order_diffusion_type = low_order_diffusion_type;
      break;
    case HighOrderScheme::entropy_diff:
      AssertThrow(false, ExcNotImplemented());
      entropy_viscosity_type = ViscosityType::none;
      high_order_viscosity_type = ViscosityType::none;
      high_order_diffusion_type = DiffusionType::entropy;
      break;
    default:
      AssertThrow(false, ExcNotImplemented());
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
  if (are_star_states)
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
  if (problem_base_parameters->boundary_conditions_type == "dirichlet")
  {
    for (unsigned int component = 0; component < n_components; ++component)
    {
      // specify to impose Dirichlet BC for all components of solution
      std::vector<bool> component_mask(n_components, true);

      problem_base_parameters->dirichlet_function->set_time(0.0);
      VectorTools::interpolate_boundary_values(
        dof_handler,
        0, // boundary ID for Dirichlet boundary
        *(problem_base_parameters->dirichlet_function),
        constraints,
        component_mask);
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
  system_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
  consistent_mass_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
  lumped_mass_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
  low_order_diffusion_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
  high_order_diffusion_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
  inviscid_ss_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
  low_order_ss_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, dsp_unconstrained, mpi_communicator);
#else
  // reinitialize matrices with sparsity pattern
  system_matrix.reinit(unconstrained_sparsity_pattern);
  consistent_mass_matrix.reinit(unconstrained_sparsity_pattern);
  lumped_mass_matrix.reinit(unconstrained_sparsity_pattern);
  low_order_diffusion_matrix.reinit(unconstrained_sparsity_pattern);
  high_order_diffusion_matrix.reinit(unconstrained_sparsity_pattern);
  inviscid_ss_matrix.reinit(unconstrained_sparsity_pattern);
  low_order_ss_matrix.reinit(unconstrained_sparsity_pattern);
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
  tmp_vector.reinit(locally_owned_dofs, mpi_communicator);
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
  tmp_vector.reinit(n_dofs);
  estimated_error_per_cell.reinit(triangulation.n_active_cells());
#endif

  // reinitialize objects
  if (are_star_states)
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
    // Lax viscosity
    case ViscosityType::lax:
    {
      // ensure diffusion type is compatible
      Assert(low_order_diffusion_type == DiffusionType::none ||
               low_order_diffusion_type == DiffusionType::laplacian,
             ExcInvalidDiffusionType());

      // create viscosity multiplier
      auto viscosity_multiplier = create_viscosity_multiplier();

      // create low-order viscosity
      low_order_viscosity =
        std::make_shared<LowOrderViscosity<dim>>(parameters.lax_viscosity_coef,
                                                 max_flux_speed_cell,
                                                 fe,
                                                 dof_handler,
                                                 cell_quadrature,
                                                 viscosity_multiplier);

      break;
    }
    // DMP viscosity
    case ViscosityType::DMP:
    {
      // ensure diffusion type is compatible
      Assert(low_order_diffusion_type == DiffusionType::none ||
               low_order_diffusion_type == DiffusionType::graphtheoretic,
             ExcInvalidDiffusionType());

      // create low-order viscosity
      low_order_viscosity = std::make_shared<DMPLowOrderViscosity<dim>>(
        dof_handler, inviscid_ss_matrix, dofs_per_cell);

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
      AssertThrow(false, ExcNotImplemented());
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
      // create entropy
      auto entropy = create_entropy();

      // set flag for if viscosity is to be used in a Laplacian
      // (vs. graph-theoretic) diffusion term
      const bool use_in_laplacian_term =
        high_order_diffusion_type == DiffusionType::laplacian;

      // create entropy viscosity
      entropy_viscosity =
        std::make_shared<EntropyViscosity<dim>>(parameters,
                                                entropy,
                                                fe,
                                                dof_handler,
                                                cell_quadrature,
                                                face_quadrature,
                                                use_in_laplacian_term);
      break;
    }
    default:
    {
      AssertThrow(false, ExcNotImplemented());
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
      AssertThrow(false, ExcNotImplemented());
      break;
    }
  }

  // create low-order and high-order artificial diffusion
  low_order_diffusion = create_artificial_diffusion(low_order_diffusion_type);
  high_order_diffusion = create_artificial_diffusion(high_order_diffusion_type);

  // get list of Dirichlet DoF indices
  get_dirichlet_values(dirichlet_values);
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
      AssertThrow(false, ExcNotImplemented());
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
    minimum_cell_diameter = std::min(minimum_cell_diameter, cell->diameter());
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
  linear_solver = std::make_shared<LinearSolver<dim>>(
    parameters,
    constraints,
    dof_handler,
    problem_base_parameters->dirichlet_function,
    n_components);

  // create SSPRK time integrator
  SSPRKTimeIntegrator<dim> ssprk(parameters.ssprk_discretization,
                                 n_dofs,
                                 *linear_solver,
                                 unconstrained_sparsity_pattern);

  // create FCT if applicable
  std::shared_ptr<ExplicitEulerFCT<dim>> fct;
  if (parameters.scheme == Scheme::fct)
    fct = create_fct();

  // if using DMP low-order viscosity, compute inviscid steady-state matrix
  bool inviscid_ss_matrix_computed = false;
  if (need_to_compute_inviscid_ss_matrix)
  {
    compute_inviscid_ss_matrix(old_solution, inviscid_ss_matrix);
    inviscid_ss_matrix_computed = true;
  }

  // if using DMP CFL condition, compute DMP CFL time step size
  double dt_dmp_cfl = 0.0;
  if (parameters.time_step_size_option == TimeStepSizeOption::cfl_dmp)
  {
    // assert that conservation law is linear; otherwise, the low-order
    // steady-state matrix would need to be computed before each time step
    AssertThrow(is_linear, ExcNotImplemented());

    // compute the inviscid steady-state matrix if it is not yet computed
    if (~inviscid_ss_matrix_computed)
      compute_inviscid_ss_matrix(old_solution, inviscid_ss_matrix);

    // compute the low-order diffusion matrix
    low_order_viscosity->update(old_solution, old_solution, 0.0, 0);
    low_order_diffusion->compute_diffusion_matrix(
      old_solution, low_order_viscosity, low_order_diffusion_matrix);

    // compute the low-order steady-state matrix
    low_order_ss_matrix.copy_from(inviscid_ss_matrix);
    low_order_ss_matrix.add(1.0, low_order_diffusion_matrix);

    // compute DMP CFL time step size
    dt_dmp_cfl = compute_dt_from_dmp_cfl_condition();
  }

  // compute end time
  double end_time;
  if (parameters.use_default_end_time &&
      problem_base_parameters->has_default_end_time)
    end_time = problem_base_parameters->default_end_time;
  else
    end_time = parameters.end_time;

  // initialize old time, time index, and transient flag
  old_time = 0.0;
  double old_dt = 1.0; // arbitrary value; shouldn't have any effect
  unsigned int n = 1;
  bool in_transient = true;

  // start transient
  while (in_transient)
  {
    // Assert that if the inviscid ss matrix needs to be computed, that it is
    // not nonlinear; otherwise it would need to be computed not just at each
    // time step, but at each stage.
    Assert(!need_to_compute_inviscid_ss_matrix || is_linear, ExcNotImplemented());

    // update max speed for use in CFL computation
    update_flux_speeds();

    // compute dt
    double dt;
    switch (parameters.time_step_size_option)
    {
      case TimeStepSizeOption::constant:
        dt = parameters.time_step_size;
        break;
      case TimeStepSizeOption::cfl:
        dt = compute_dt_from_cfl_condition();
        break;
      case TimeStepSizeOption::cfl_dmp:
        dt = dt_dmp_cfl;
        break;
      case TimeStepSizeOption::cfl_di:
        dt = compute_dt_from_di_cfl_condition();
        break;
      default:
        AssertThrow(false, ExcNotImplemented());
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

    // preserve default outstream flags
    std::ios::fmtflags default_flags(std::cout.flags());

    // print CFL in red if it violates CFL condition
    if (cfl_is_violated)
      cout1 << std::fixed << std::setprecision(5) << "  time step " << n
            << ": t = \x1b[1;34m" << new_time << "\x1b[0m, CFL = \x1b[1;31m"
            << cfl << "\x1b[0m" << std::endl;
    else
      cout1 << std::fixed << std::setprecision(5) << "  time step " << n
            << ": t = \x1b[1;34m" << new_time << "\x1b[0m, CFL = " << cfl
            << std::endl;

    // restore default outstream flags
    std::cout.flags(default_flags);

    // initialize SSPRK time step
    ssprk.initialize_time_step(old_solution, dt);

    // loop over SSPRK stages to compute new solution
    for (unsigned int i = 0; i < ssprk.n_stages; ++i)
    {
      // print stage index
      cout1 << "    stage " << i + 1 << " of " << ssprk.n_stages << std::endl;

      // get stage time
      double t_stage = ssprk.get_stage_time();

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

      // compute steady-state rhs vector
      compute_ss_rhs(t_stage, ss_rhs);

      // compute steady-state flux vector
      compute_ss_flux(dt, old_stage_solution, ss_flux);

      // advance by an SSPRK step
      switch (parameters.scheme)
      {
        case Scheme::low:
          low_order_viscosity->update(
            old_stage_solution, older_solution, old_stage_dt, n);
          low_order_diffusion->compute_diffusion_matrix(
            old_stage_solution, low_order_viscosity, low_order_diffusion_matrix);
          ssprk.step(lumped_mass_matrix,
                     ss_flux,
                     low_order_diffusion_matrix,
                     ss_rhs,
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
                     high_order_diffusion_matrix,
                     ss_rhs,
                     true);
          break;
        case Scheme::fct:
          perform_fct_ssprk_step(dt, old_stage_dt, n, fct, t_stage, ssprk);
          break;
        default:
          AssertThrow(false, ExcNotImplemented());
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
      if (fct && parameters.output_transient_fct_bounds)
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

  // FCT output
  if (fct)
  {
    // output final FCT bounds if specified
    if (parameters.output_final_fct_bounds)
      fct->output_bounds(postprocessor);

    // output limiter matrix if specified
    if (parameters.output_limiter_matrix)
      fct->output_limiter_matrix();
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

/**
 * \brief Computes time step size using the DMP CFL condition.
 *
 * The DMP CFL condition for the explicit Euler method is the following:
 * \f[
 *   \Delta t^n \leq \frac{M^L_{i,i}}{A^{L,n}_{i,i}} \quad \forall i \,,
 * \f]
 * where \f$\mathbf{M}^L\f$ is the lumped mass matrix,
 * and \f$\mathbf{A}^{L,n}\f$ is the low-order steady-state matrix evaluated
 * at time \f$n\f$. This function computes the time step size as
 * \f[
 *   \Delta t^n = \nu\min\limits_i\frac{M^L_{i,i}}{A^{L,n}_{i,i}} \,,
 * \f]
 * where \f$\nu\f$ is the user-specified CFL number.
 * Dirichlet nodes are excluded from the \f$\min\f$ because the time step
 * restriction does not apply and the condition \f$A^{L,n} > 0\f$ is not
 * necessarily satisfied.
 *
 * \pre This function assumes that the lumped mass matrix and the low-order
 * steady-state matrix have been computed. If the low-order steady-state
 * matrix depends on the solution or explicitly on time, this function
 * assumes it has computed for this time step.
 */
template <int dim>
double ConservationLaw<dim>::compute_dt_from_dmp_cfl_condition()
{
  // iterators for determining if an index is in the list of Dirichlet indices
  std::map<unsigned int, double>::iterator it_end = dirichlet_values.end();

  // initialize CFL time step size to arbitrary large number before min()
  double dt = 1.0e15;
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    // if not a Dirichlet node
    if (dirichlet_values.find(i) == it_end)
    {
      const double ALii = low_order_ss_matrix(i, i);
      Assert(ALii > 0.0, ExcNegativeDiagonal(i, ALii));
      dt = std::min(dt, lumped_mass_matrix(i, i) / ALii);
    }
  }

  // multiply by user-specified CFL number
  dt *= parameters.cfl;

  return dt;
}

/**
 * \brief Computes time step size using the DI CFL condition.
 */
template <int dim>
double ConservationLaw<dim>::compute_dt_from_di_cfl_condition()
{
  AssertThrow(false, ExcNotImplemented());

  return 0.0;
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
    viscosities.push_back(entropy_viscosity);
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
 * \param[in] t_stage  stage time
 * \param[in,out] ssprk SSPRK time integrator
 */
template <int dim>
void ConservationLaw<dim>::perform_fct_ssprk_step(
  const double & dt,
  const double & old_stage_dt,
  const unsigned int & n,
  const std::shared_ptr<ExplicitEulerFCT<dim>> & fct,
  const double & t_stage,
  SSPRKTimeIntegrator<dim> & ssprk)
{
  // update star states
  if (are_star_states)
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
    consistent_mass_matrix, ss_flux, high_order_diffusion_matrix, ss_rhs, false);

  // get high-order solution
  ssprk.get_intermediate_solution(new_solution);

  // perform FCT
  // form rhs: system_rhs = M*u_old + dt*(ss_rhs - ss_flux - D*u_old + f)
  system_rhs = 0;
  lumped_mass_matrix.vmult(tmp_vector, old_stage_solution);
  system_rhs.add(1.0, tmp_vector);
  system_rhs.add(dt, ss_rhs);
  system_rhs.add(-dt, ss_flux);
  low_order_diffusion_matrix.vmult(tmp_vector, old_stage_solution);
  system_rhs.add(-dt, tmp_vector);

  // compute antidiffusion vector
  Vector<double> & antidiffusion_vector = tmp_vector;
  fct->compute_antidiffusion_vector(new_solution,
                                    old_stage_solution,
                                    dt,
                                    ss_flux,
                                    ss_reaction,
                                    low_order_diffusion_matrix,
                                    high_order_diffusion_matrix,
                                    ss_rhs,
                                    t_stage,
                                    antidiffusion_vector);

  // add antidiffusion vector
  system_rhs.add(dt, antidiffusion_vector);

  // solve the linear system M*u_new = system_rhs
  system_matrix.copy_from(lumped_mass_matrix);
  linear_solver->solve_with_dirichlet(system_matrix, new_solution, system_rhs);

  // set stage solution to be FCT solution for this stage
  ssprk.set_intermediate_solution(new_solution);

  // check FCT bounds if specified
  if (parameters.check_fct_bounds)
    fct->check_bounds(new_solution);

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
void ConservationLaw<dim>::get_dirichlet_values(
  std::map<unsigned int, double> & values)
{
  if (problem_base_parameters->boundary_conditions_type == "dirichlet")
  {
    // loop over components
    for (unsigned int component = 0; component < n_components; ++component)
    {
      // mask other components
      std::vector<bool> component_mask(n_components, false);
      component_mask[component] = true;

      // fill boundary_values with boundary values
      VectorTools::interpolate_boundary_values(
        dof_handler,
        0, // boundary ID for Dirichlet boundary
        *(problem_base_parameters->dirichlet_function),
        values,
        component_mask);
    }
  }
}

/**
 * \brief Computes the inviscid steady state matrix \f$\mathbf{A}\f$.
 *
 * \param[in] solution  solution vector
 * \param[out] matrix   inviscid steady-state matrix
 */
template <int dim>
void ConservationLaw<dim>::compute_inviscid_ss_matrix(const Vector<double> &,
                                                      SparseMatrix<double> &)
{
  // throw exception if this base version is called; if a derived class needs
  // this function (i.e., it is a scalar conservation law), it should override
  // this function
  AssertThrow(false, ExcInvalidState());
}

/**
 * \brief Computes the steady-state reaction vector.
 *
 * This function computes the steady-state reaction vector, which has the
 * entries
 * \f[
 *   \sigma_i = \int\limits_{S_i}\varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *     = \sum\limits_j\int\limits_{S_{i,j}}
 *     \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \,,
 * \f]
 * where \f$\sigma(\mathbf{x})\f$ is the reaction coefficient of the
 * conservation law equation when it is put in the form
 * \f[
 *   \frac{\partial\mathbf{u}}{\partial t}
 *   + \nabla \cdot \mathbf{F}(\mathbf{u}) + \sigma(\mathbf{x})\mathbf{u}
 *   = \mathbf{0} \,.
 * \f]
 *
 * \param[out] ss_reaction steady-state reaction vector
 */
template <int dim>
void ConservationLaw<dim>::compute_ss_reaction(Vector<double> & ss_reaction)
{
  ss_reaction = 0;
}

/**
 * \brief Creates a star state object and returns the pointer.
 *
 * This function throws an exception in this base class.
 *
 * \return star state object pointer
 */
template <int dim>
std::shared_ptr<StarState<dim>> ConservationLaw<dim>::create_star_state() const
{
  // throw exception if not overridden by derived class
  AssertThrow(false, ExcNotImplemented());

  // return null pointer to avoid compiler warning about not returning anything
  return nullptr;
}
