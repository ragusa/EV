/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 */
template <int dim>
TransportProblem<dim>::TransportProblem(
  const TransportRunParameters & run_parameters_)
  : cout1(std::cout, run_parameters_.verbosity_level >= 1),
    cout2(std::cout, run_parameters_.verbosity_level >= 2),
    run_parameters(run_parameters_),
    is_time_dependent(!(run_parameters_.temporal_discretization ==
                        TemporalDiscretizationClassification::ss)),
    problem_parameters(run_parameters_.problem_name, !is_time_dependent),
    timer(cout2, TimerOutput::summary, TimerOutput::wall_times)
{
}

/**
 * \brief Initializes system.
 */
template <int dim>
void TransportProblem<dim>::initialize_system()
{
  // timer
  TimerOutput::Scope t_initialize(timer, "initialize");

  // get problem parameters
  get_problem_parameters();

  // initially refine mesh
  triangulation.refine_global(run_parameters.initial_refinement_level);
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
    source_path + "/problems/transport/" + run_parameters.problem_name;

  // get and process the problem parameters
  const FESystem<dim> fe(FE_Q<dim>(run_parameters.degree), 1);
  const QGauss<dim - 1> face_quadrature(run_parameters.n_quadrature_points);
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
  initialize_system();

  // determine end time; either use user-specified value or default
  double end_time;
  if (run_parameters.use_default_end_time &&
      problem_parameters.has_default_end_time)
    end_time = problem_parameters.default_end_time;
  else
    end_time = run_parameters.end_time;

  // solution component names and interpretations
  std::vector<std::string> component_names(1, "angular_flux");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretations(1,
                              DataComponentInterpretation::component_is_scalar);

  // create post-processor object
  PostProcessor<dim> postprocessor(run_parameters,
                                   1,
                                   end_time,
                                   problem_parameters.has_exact_solution,
                                   problem_parameters.exact_solution_function,
                                   run_parameters.problem_name,
                                   component_names,
                                   component_interpretations,
                                   triangulation);

  // create refinement handler object
  RefinementHandler<dim> refinement_handler(run_parameters, triangulation);

  // loop over refinement cycles
  for (unsigned int cycle = 0; cycle < run_parameters.n_refinement_cycles;
       ++cycle)
  {
    // set cycle for post-processor
    postprocessor.set_cycle(cycle);

    // refine
    refinement_handler.refine(cycle);

    // print information
    cout1 << std::endl << "Cycle " << cycle << "  (n_cells = ";
    cout1 << triangulation.n_active_cells() << ")" << std::endl;

    if (is_time_dependent)
    { // run transient problem
      // timer
      TimerOutput::Scope t_solve(timer, "solve");

      // get nominal time step size for cycle
      const double nominal_dt = refinement_handler.get_nominal_time_step_size();

      // create and run appropriate executioner
      switch (run_parameters.temporal_discretization)
      {
        case TemporalDiscretizationClassification::ssprk:
        {
          TransportSSPRKExecutioner<dim> executioner(run_parameters,
                                                     problem_parameters,
                                                     triangulation,
                                                     postprocessor,
                                                     nominal_dt);
          executioner.run();
          break;
        }
        case TemporalDiscretizationClassification::theta:
        {
          throw ExcNotImplemented();
          break;
        }
        default:
        {
          throw ExcNotImplemented();
          break;
        }
      }
    }
    else
    { // run steady-state problem
      // timer
      TimerOutput::Scope t_solve(timer, "solve");

      // create and run steady-state executioner
      TransportSteadyStateExecutioner<dim> executioner(
        run_parameters, problem_parameters, triangulation, postprocessor);
      executioner.run();
    }
  } // end refinement cycle loop
}
