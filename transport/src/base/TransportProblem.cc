/**
 * \brief Constructor.
 *
 * \param[in] parameters_  input parameters
 */
template <int dim>
TransportProblem<dim>::TransportProblem(
  const TransportParameters<dim> & parameters_)
  : parameters(parameters_),
    is_time_dependent(
      !(parameters.temporal_discretization == TemporalDiscretization::ss)),
    problem_parameters(parameters.problem_name, !is_time_dependent),
    timer(std::cout, TimerOutput::summary, TimerOutput::wall_times)
{
}

/**
 * \brief Initializes system.
 */
template <int dim>
void TransportProblem<dim>::initializeSystem()
{
  // timer
  TimerOutput::Scope t_initialize(timer, "initialize");

  // get problem parameters
  get_problem_parameters();

  // initially refine mesh
  triangulation.refine_global(parameters.initial_refinement_level);
}

/**
 * \brief Gets problem parameters.
 */
template <int dim>
void TransportProblem<dim>::get_problem_parameters()
{
  // get path of source directory from #define created by CMake
  std::stringstream source_path_ss;
  source_path_ss << SOURCE_PATH;
  std::string source_path;
  source_path_ss >> source_path;

  // create problem parameters file name
  std::string problem_parameters_file =
    source_path + "/problems/transport/" + parameters.problem_name;

  // get and process the problem parameters
  const FESystem<dim> fe(FE_Q<dim>(parameters.degree), 1);
  const QGauss<dim - 1> face_quadrature(parameters.n_quadrature_points);
  problem_parameters.get_and_process_parameters(
    problem_parameters_file, triangulation, fe, face_quadrature);
}

/**
 * \brief Runs the problem.
 */
template <int dim>
void TransportProblem<dim>::run()
{
  // initialize problem
  initializeSystem();

  // create post-processor object
  PostProcessor<dim> postprocessor(parameters,
                                   problem_parameters.has_exact_solution,
                                   problem_parameters.exact_solution_function);

  // create refinement handler object
  RefinementHandler<dim> refinement_handler(parameters, triangulation);

  // loop over refinement cycles
  for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; ++cycle)
  {
    // set cycle for post-processor
    postprocessor.setCycle(cycle);

    // refine
    refinement_handler.refine(cycle);

    // print information
    std::cout << std::endl << "Cycle " << cycle << "  (n_cells = ";
    std::cout << triangulation.n_active_cells() << ")" << std::endl;

    if (is_time_dependent)
    { // run transient problem
      // timer
      TimerOutput::Scope t_solve(timer, "solve");

      // get nominal time step size for cycle
      const double nominal_dt = refinement_handler.get_nominal_time_step_size();

      // create and run transient executioner
      TransientExecutioner<dim> executioner(
        parameters,
        triangulation,
        problem_parameters.transport_direction,
        problem_parameters.transport_speed,
        problem_parameters.cross_section_function,
        problem_parameters.source_function,
        *(problem_parameters.dirichlet_function),
        problem_parameters.initial_conditions_function,
        problem_parameters.domain_volume,
        postprocessor,
        problem_parameters.source_is_time_dependent,
        nominal_dt);
      executioner.run();
    }
    else
    { // run steady-state problem
      // timer
      TimerOutput::Scope t_solve(timer, "solve");

      // create and run steady-state executioner
      SteadyStateExecutioner<dim> executioner(
        parameters,
        triangulation,
        problem_parameters.transport_direction,
        problem_parameters.transport_speed,
        problem_parameters.cross_section_function,
        problem_parameters.source_function,
        *(problem_parameters.dirichlet_function),
        problem_parameters.domain_volume,
        postprocessor);
      executioner.run();
    }
  } // end refinement cycle loop
}

