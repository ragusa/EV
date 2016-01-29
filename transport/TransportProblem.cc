/**
 * Constructor.
 *
 * @param[in] parameters input parameters
 */
template<int dim>
TransportProblem<dim>::TransportProblem(
  const TransportParameters<dim> &parameters) :
    parameters(parameters),
    is_time_dependent(
      !(parameters.time_discretization_option == TransportParameters<dim>::SS)),
    transport_direction(),
    exact_solution_option(ExactSolutionOption::none),
    has_exact_solution(false),
    timer(std::cout, TimerOutput::summary, TimerOutput::wall_times)
{
}

/**
 * Destructor.
 */
template<int dim>
TransportProblem<dim>::~TransportProblem()
{
}

/**
 * Initializes system.
 */
template<int dim>
void TransportProblem<dim>::initializeSystem()
{
  // timer
  TimerOutput::Scope t_initialize(timer, "initialize");

  // process problem ID
  processProblemID();

  // determine function variables based on dimension;
  // this is to be used in initialization of function parser objects
  std::string variables;
  if (dim == 1)
    variables = "x";
  else if (dim == 2)
    variables = "x,y";
  else if (dim == 3)
    variables = "x,y,z";
  else
    Assert(false, ExcInvalidState());

  // add time variable if not steady state
  if (is_time_dependent)
    variables += ",t";

  // create function parser constants
  function_parser_constants["pi"] = numbers::PI;
  function_parser_constants["x_min"] = x_min;
  function_parser_constants["x_mid"] = x_min + 0.5 * (x_max - x_min);
  function_parser_constants["x_max"] = x_max;

  // set flag for problem having exact solution
  has_exact_solution = exact_solution_option != ExactSolutionOption::none;

  // initialize exact solution function
  if (has_exact_solution)
  {
    if (exact_solution_option == ExactSolutionOption::parser)
    { // function parser

      // create and initialize function parser
      std::shared_ptr<FunctionParser<dim> > exact_solution_function_derived =
        std::make_shared<FunctionParser<dim> >();
      exact_solution_function_derived->initialize(variables,
        exact_solution_string, function_parser_constants, is_time_dependent);

      // point base class shared pointer to derived class function object
      exact_solution_function = exact_solution_function_derived;

    }
    else if (exact_solution_option != ExactSolutionOption::multi_region)
    { // multi-region
      ExcInvalidState();
    }
  }

  // initialize initial conditions function
  if (is_time_dependent)
    initial_conditions.initialize(variables, initial_conditions_string,
      function_parser_constants, is_time_dependent);

  // initialize source function
  source_function.initialize(variables, source_string, function_parser_constants,
    is_time_dependent);

  // initialize cross section function
  cross_section_function.initialize(variables, cross_section_string,
    function_parser_constants, is_time_dependent);

  // initialize Dirichlet boundary value function
  incoming_function.initialize(variables, incoming_string,
    function_parser_constants, is_time_dependent);

  // create grid for initial refinement level
  GridGenerator::hyper_cube(triangulation, x_min, x_max);
  triangulation.refine_global(parameters.initial_refinement_level);

  // compute domain volume - assume hypercube
  domain_volume = std::pow(x_max - x_min, dim);
}

/**
 * Processes problem ID.
 */
template<int dim>
void TransportProblem<dim>::processProblemID()
{
  switch (parameters.problem_id)
  {
    case 1 :
    { // pure absorber
      Assert(dim < 3, ExcNotImplemented());

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      incoming_string = "1";
      function_parser_constants["incoming"] = 1.0;

      cross_section_string = "1.0";
      function_parser_constants["sigma"] = 1.0;

      source_time_dependent = false;
      source_string = "0";
      function_parser_constants["source"] = 0.0;

      exact_solution_option = ExactSolutionOption::parser;

      if (!is_time_dependent)
        exact_solution_string =
          "source/sigma + (incoming - source/sigma)*exp(-sigma*(x-x_min))";
      else
        exact_solution_string =
          "if(x<=t,source/sigma + (incoming - source/sigma)*exp(-sigma*(x-x_min)),0)";

      initial_conditions_string = "0";

      break;
    }
    case 2 :
    { // void-to-absorber
      Assert(dim < 3, ExcNotImplemented());

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      incoming_string = "1";
      function_parser_constants["incoming"] = 1.0;

      if (dim == 1)      // 1-D
        cross_section_string = "if(x<x_mid, 0, sigma)";
      else if (dim == 2) // 2-D
        cross_section_string = "if(x>=x_mid, if(y>=x_mid, sigma, 0), 0)";
      else
        // 3-D
        cross_section_string = "if(x>=x_mid, if(y>=x_mid, if(z>=x_mid,"
          "sigma, 0), 0), 0)";
      function_parser_constants["sigma"] = 100.0;

      source_time_dependent = false;
      source_string = "0";
      function_parser_constants["source"] = 0.0;

      exact_solution_option = ExactSolutionOption::parser;

      if (!is_time_dependent)
      { // steady-state
        if (dim == 1)      // 1-D
          exact_solution_string = "if(x>=x_mid,"
            "incoming*exp(-sigma*(x-x_mid)), incoming)";
        else if (dim == 2) // 2-D
          exact_solution_string = "if(x>=x_mid, if(y>=y_mid,"
            "incoming*exp(-sigma*(x-x_mid)), incoming), incoming)";
        else
          // 3-D
          exact_solution_string = "if(x>=x_mid, if(y>=y_mid, if(z>=x_mid,"
            "incoming*exp(-sigma*(x-x_mid)), incoming), incoming),"
            "incoming)";
      }
      else
      { // transient
        if (dim == 1)      // 1-D
          exact_solution_string = "if(x-t<x_min, if(x>=x_mid,"
            "incoming*exp(-sigma*(x-x_mid)), incoming), 0)";
        else if (dim == 2) // 2-D
          exact_solution_string = "if(x-t<x_min, if(x>=x_mid,"
            "if(y>=x_mid, incoming*exp(-sigma*(x-x_mid)), incoming),"
            "incoming), 0)";
        else
          // 3-D
          exact_solution_string = "if(x-t<x_min, if(x>=x_mid,"
            "if(y>=x_mid, if(z>=x_mid, incoming*exp(-sigma*(x-x_mid)),"
            "incoming), incoming), incoming), 0)";
      }
      initial_conditions_string = "0";
      break;
    }
    case 3 :
    { // void
      Assert(dim == 1, ExcNotImplemented());

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      cross_section_string = "0";
      source_string = "0";
      source_time_dependent = false;

      incoming_string = "1";
      function_parser_constants["incoming"] = 1.0;

      exact_solution_option = ExactSolutionOption::parser;

      if (!is_time_dependent)
        exact_solution_string = "1.0";
      else
        exact_solution_string = "if((x-t)<0,incoming,0)";
      initial_conditions_string = "if(x<0,incoming,0)";
      break;
    }
    case 5 :
    { // MMS-1
      Assert(dim == 1, ExcNotImplemented()); // assume 1-D

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      incoming_string = "0";

      cross_section_string = "1.0";

      exact_solution_option = ExactSolutionOption::parser;

      // use different MS for steady-state vs. transient problem
      if (is_time_dependent)
      {
        exact_solution_string = "t*sin(pi*x)"; // assume omega_x = 1 and c = 1
        source_string = "sin(pi*x) + pi*t*cos(pi*x) + t*sin(pi*x)";
        source_time_dependent = true;
        initial_conditions_string = exact_solution_string;
      }
      else
      {
        exact_solution_string = "sin(pi*x)"; // assume omega_x = 1 and c = 1
        source_string = "pi*cos(pi*x) + sin(pi*x)";
        source_time_dependent = false;
      }

      break;
    }
    case 6 :
    { // MMS-2
      Assert(dim == 1, ExcNotImplemented()); // assume 1-D
      Assert(is_time_dependent, ExcNotImplemented()); // assume not steady-state

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      incoming_string = "0";

      cross_section_string = "1";

      exact_solution_option = ExactSolutionOption::parser;
      exact_solution_string = "x*exp(-t)"; // assume omega_x = 1 and c = 1

      source_string = "-x*exp(-t) + exp(-t) + x*exp(-t)";
      source_time_dependent = true;

      initial_conditions_string = exact_solution_string;

      break;
    }
    case 7 :
    { // MMS-3
      Assert(dim == 1, ExcNotImplemented()); // assume 1-D
      Assert(is_time_dependent, ExcNotImplemented()); // assume not steady-state

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      incoming_string = "0";

      cross_section_string = "1";

      exact_solution_option = ExactSolutionOption::parser;
      exact_solution_string = "exp(-t)*sin(pi*x)"; // assume omega_x = 1 and c = 1

      source_string =
        "-exp(-t)*sin(pi*x) + pi*exp(-t)*cos(pi*x) + exp(-t)*sin(pi*x)";
      source_time_dependent = true;

      initial_conditions_string = exact_solution_string;

      break;
    }
    case 8 :
    { // source in left half
      Assert(dim == 1, ExcNotImplemented()); // assume 1-D
      Assert(is_time_dependent, ExcNotImplemented()); // assume not steady-state

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      function_parser_constants["speed"] = 1.0;

      incoming_string = "0";
      function_parser_constants["incoming"] = 0.0;

      cross_section_string = "100";
      function_parser_constants["sigma"] = 100.0;

      source_time_dependent = false;
      source_string = "if (x<x_mid,10,0)";
      function_parser_constants["source"] = 10.0;

      exact_solution_option = ExactSolutionOption::parser;
      exact_solution_string = (std::string) "source/sigma*(1-exp(-sigma*"
        + "max(0,min(x,x_mid)-max(x-speed*t,0))))"
        + "*exp(-sigma*max(0,min(x,x_max)-max(x-speed*t,x_mid)))";

      initial_conditions_string = "0";

      break;
    }
    case 9 :
    { // MMS-4
      Assert(dim == 1, ExcNotImplemented()); // assume 1-D
      Assert(is_time_dependent, ExcNotImplemented()); // assume not steady-state

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      incoming_string = "0";

      cross_section_string = "1";

      exact_solution_option = ExactSolutionOption::parser;
      exact_solution_string = "x*t"; // assume omega_x = 1 and c = 1

      source_string = "x + t + x*t";
      source_time_dependent = true;

      initial_conditions_string = exact_solution_string;

      break;
    }
    case 10 :
    { // MMS-5
      Assert(dim == 1, ExcNotImplemented()); // assume 1-D

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      cross_section_string = "1";

      exact_solution_option = ExactSolutionOption::parser;

      if (is_time_dependent)
      {
        incoming_string = "t";
        exact_solution_string = "t"; // assume omega_x = 1 and c = 1
        source_string = "1 + t";
        source_time_dependent = true;
        initial_conditions_string = exact_solution_string;
      }
      else
      {
        incoming_string = "1";
        exact_solution_string = "1"; // assume omega_x = 1 and c = 1
        source_string = "1";
        source_time_dependent = false;
      }

      break;
    }
    case 11 :
    { // skew void-to-absorber

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0 / sqrt(2.0);
      if (dim >= 2)
        transport_direction[1] = 1.0 / sqrt(3.0);
      if (dim >= 3)
        transport_direction[2] = 1.0 / sqrt(6.0);

      incoming_string = "1";
      function_parser_constants["incoming"] = 1.0;

      if (dim == 1)      // 1-D
        cross_section_string = "if(x<x_mid, 0, sigma)";
      else if (dim == 2) // 2-D
        cross_section_string = "if(x>=x_mid, if(y>=x_mid, sigma, 0), 0)";
      else
        // 3-D
        cross_section_string = "if(x>=x_mid, if(y>=x_mid, if(z>=x_mid,"
          "sigma, 0), 0), 0)";
      function_parser_constants["sigma"] = 10.0;

      source_time_dependent = false;
      source_string = "0";
      function_parser_constants["source"] = 0.0;

      // for now, assume no steady-state, but this would be easy to implement
      // by using the existing transient exact solution class and just using
      // t=large
      Assert(is_time_dependent, ExcNotImplemented());

      // create exact solution function constructor arguments
      const std::vector<double> interface_positions = { 0.5 };
      const std::vector<double> region_sources = { 0.0, 0.0 };
      const std::vector<double> region_sigmas = { 0.0, 10.0 };
      const std::vector<double> direction(
        { 1.0 / sqrt(2.0), 1.0 / sqrt(3.0), 1.0 / sqrt(6.0) });
      const double incoming = 1.0;
      exact_solution_option = ExactSolutionOption::multi_region;

      // create MultiRegionExactSolution object
      std::shared_ptr<MultiRegionExactSolution<dim> > exact_solution_function_derived =
        std::make_shared<MultiRegionExactSolution<dim> >(interface_positions,
          region_sources, region_sigmas, direction, incoming);
      // point base class shared pointer to derived class function object
      exact_solution_function = exact_solution_function_derived;

      initial_conditions_string = "0";

      break;
    }
    case 12 :
    { // 3-region

      // for now, assume 1-D
      Assert(dim == 1, ExcNotImplemented());

      x_min = 0.0;
      x_max = 1.0;

      transport_direction[0] = 1.0;

      incoming_string = "1";
      function_parser_constants["incoming"] = 1.0;

      cross_section_string = "if(x<=x1, sigma0, if(x<=x2, sigma1, sigma2))";
      function_parser_constants["sigma0"] = 1.0;
      function_parser_constants["sigma1"] = 40.0;
      function_parser_constants["sigma2"] = 20.0;
      function_parser_constants["x1"] = 0.3;
      function_parser_constants["x2"] = 0.6;

      source_time_dependent = false;
      source_string = "if(x<=x1, q0, if(x<=x2, q1, q2))";
      function_parser_constants["q0"] = 1.0;
      function_parser_constants["q1"] = 5.0;
      function_parser_constants["q2"] = 20.0;

      // create exact solution function constructor arguments
      const std::vector<double> interface_positions = { 0.3, 0.6 };
      const std::vector<double> region_sources = { 1.0, 5.0, 20.0 };
      const std::vector<double> region_sigmas = { 1.0, 40.0, 20.0 };
      const std::vector<double> direction( { 1.0, 0.0, 0.0 });
      const double incoming = 1.0;
      exact_solution_option = ExactSolutionOption::multi_region;

      // create MultiRegionExactSolution object
      std::shared_ptr<MultiRegionExactSolution<dim> > exact_solution_function_derived =
        std::make_shared<MultiRegionExactSolution<dim> >(interface_positions,
          region_sources, region_sigmas, direction, incoming);
      // point base class shared pointer to derived class function object
      exact_solution_function = exact_solution_function_derived;

      // for now, assume no steady-state, but this would be easy to implement
      // by using the existing transient exact solution class and just using
      // t=large
      Assert(is_time_dependent, ExcNotImplemented());

      initial_conditions_string = "0";

      break;
    }
    case 13 :
    { // source-in-void to absorber

      // create exact solution function constructor arguments
      const std::vector<double> interface_positions = { 0.5 };
      const std::vector<double> region_sources = { 1.0, 0.0 };
      const std::vector<double> region_sigmas = { 0.0, 10.0 };
      const std::vector<double> direction( { 1.0, 0.0, 0.0 });
      const double incoming = 0.0;

      x_min = 0.0;
      x_max = 1.0;

      for (unsigned int d = 0; d < dim; ++d)
        transport_direction[d] = direction[d];

      incoming_string = "incoming";
      function_parser_constants["incoming"] = incoming;

      function_parser_constants["x1"] = interface_positions[0];

      function_parser_constants["sigma0"] = region_sigmas[0];
      function_parser_constants["sigma1"] = region_sigmas[1];

      source_time_dependent = false;
      function_parser_constants["q0"] = region_sources[0];
      function_parser_constants["q1"] = region_sources[1];

      if (dim == 1)
      {
        cross_section_string = "if(x<=x1, sigma0, sigma1)";
        source_string = "if(x<=x1, q0, q1)";
      }
      else if (dim == 2)
      {
        cross_section_string = "if(x<=x1 || y<=x1, sigma0, sigma1)";
        source_string = "if(x<=x1 || y<=x1, q0, q1)";
      }
      else
      {
        Assert(false, ExcNotImplemented());
      }

      // create MultiRegionExactSolution object
      if (is_time_dependent)
      {
        exact_solution_option = ExactSolutionOption::multi_region;

        std::shared_ptr<MultiRegionExactSolution<dim> > exact_solution_function_derived =
          std::make_shared<MultiRegionExactSolution<dim> >(interface_positions,
            region_sources, region_sigmas, direction, incoming);
        // point base class shared pointer to derived class function object
        exact_solution_function = exact_solution_function_derived;
      }
      else
      {
        exact_solution_option = ExactSolutionOption::parser;

        Assert(dim == 1, ExcNotImplemented());
        exact_solution_string = "if(x<0.5,x,0.5*exp(-10*(x-0.5)))";
      }

      initial_conditions_string = "0";

      break;
    }
    case 14 :
    { // multi-region

      const std::vector<double> interface_positions = { 0.2, 0.4, 0.6, 0.8 };
      const std::vector<double> region_sources = { 0.0, 0.0, 10.0, 10.0, 5.0 };
      const std::vector<double> region_sigmas = { 0.0, 10.0, 0.0, 10.0, 5.0 };
      const std::vector<double> direction(
        { 1.0 / sqrt(2.0), 1.0 / sqrt(3.0), 1.0 / sqrt(6.0) });
      const double incoming = 1.0;

      x_min = 0.0;
      x_max = 1.0;

      for (unsigned int d = 0; d < dim; ++d)
        transport_direction[d] = direction[d];

      incoming_string = "incoming";
      function_parser_constants["incoming"] = incoming;

      function_parser_constants["x1"] = interface_positions[0];
      function_parser_constants["x2"] = interface_positions[1];
      function_parser_constants["x3"] = interface_positions[2];
      function_parser_constants["x4"] = interface_positions[3];
      function_parser_constants["sigma0"] = region_sigmas[0];
      function_parser_constants["sigma1"] = region_sigmas[1];
      function_parser_constants["sigma2"] = region_sigmas[2];
      function_parser_constants["sigma3"] = region_sigmas[3];
      function_parser_constants["sigma4"] = region_sigmas[4];
      function_parser_constants["q0"] = region_sources[0];
      function_parser_constants["q1"] = region_sources[1];
      function_parser_constants["q2"] = region_sources[2];
      function_parser_constants["q3"] = region_sources[3];
      function_parser_constants["q4"] = region_sources[4];

      cross_section_string = "if(x<=x1 || y<=x1, sigma0, if(x<=x2 || y<=x2,"
        " sigma1, if(x<=x3 || y<=x3, sigma2, if(x<=x4 || y<=x4,"
        " sigma3, sigma4))))";

      source_time_dependent = false;
      source_string = "if(x<=x1 || y<=x1, q0, if(x<=x2 || y<=x2,"
        " q1, if(x<=x3 || y<=x3, q2, if(x<=x4 || y<=x4,"
        " q3, q4))))";

      // create exact solution function constructor arguments
      exact_solution_option = ExactSolutionOption::multi_region;

      // create MultiRegionExactSolution object
      std::shared_ptr<MultiRegionExactSolution<dim> > exact_solution_function_derived =
        std::make_shared<MultiRegionExactSolution<dim> >(interface_positions,
          region_sources, region_sigmas, direction, incoming);
      // point base class shared pointer to derived class function object
      exact_solution_function = exact_solution_function_derived;

      // for now, assume no steady-state, but this would be easy to implement
      // by using the existing transient exact solution class and just using
      // t=large
      Assert(is_time_dependent, ExcNotImplemented());

      initial_conditions_string = "0";

      break;
    }
    case 15 :
    { // 2-region saturation

      Assert(dim == 1, ExcNotImplemented());

      const std::vector<double> interface_positions = { 0.5 };
      const std::vector<double> region_sources = { 1.0, 1.0e5 };
      const std::vector<double> region_sigmas = { 1.0, 1.0e5 };
      const std::vector<double> direction( { 1.0, 0.0, 0.0 });
      const double incoming = 1.0;

      x_min = 0.0;
      x_max = 1.0;

      for (unsigned int d = 0; d < dim; ++d)
        transport_direction[d] = direction[d];

      incoming_string = "incoming";
      function_parser_constants["incoming"] = incoming;

      function_parser_constants["x1"] = interface_positions[0];
      function_parser_constants["sigma0"] = region_sigmas[0];
      function_parser_constants["sigma1"] = region_sigmas[1];
      function_parser_constants["q0"] = region_sources[0];
      function_parser_constants["q1"] = region_sources[1];

      cross_section_string = "if(x<=x1, sigma0, sigma1)";

      source_time_dependent = false;
      source_string = "if(x<=x1, q0, q1)";

      // create exact solution function constructor arguments
      exact_solution_option = ExactSolutionOption::multi_region;

      // create MultiRegionExactSolution object
      std::shared_ptr<MultiRegionExactSolution<dim> > exact_solution_function_derived =
        std::make_shared<MultiRegionExactSolution<dim> >(interface_positions,
          region_sources, region_sigmas, direction, incoming);
      // point base class shared pointer to derived class function object
      exact_solution_function = exact_solution_function_derived;

      // for now, assume no steady-state, but this would be easy to implement
      // by using the existing transient exact solution class and just using
      // t=large
      Assert(is_time_dependent, ExcNotImplemented());

      initial_conditions_string = "0";

      break;
    }
    default :
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }
}

/**
 * Runs the problem.
 */
template<int dim>
void TransportProblem<dim>::run()
{
  // initialize problem
  initializeSystem();

  // create post-processor object
  PostProcessor<dim> postprocessor(parameters, has_exact_solution,
    exact_solution_function);

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
    std::cout << std::endl << "Cycle " << cycle << ':' << std::endl;
    std::cout << "   Number of active cells:       "
      << triangulation.n_active_cells() << std::endl;

    if (is_time_dependent)
    { // run transient problem
      // timer
      TimerOutput::Scope t_solve(timer, "solve");

      // create and run transient executioner
      TransientExecutioner<dim> executioner(parameters, triangulation,
        transport_direction, cross_section_function, source_function,
        incoming_function, initial_conditions, domain_volume, postprocessor,
        source_time_dependent);
      executioner.run();
    }
    else
    { // run steady-state problem
      // timer
      TimerOutput::Scope t_solve(timer, "solve");

      // create and run steady-state executioner
      SteadyStateExecutioner<dim> executioner(parameters, triangulation,
        transport_direction, cross_section_function, source_function,
        incoming_function, domain_volume, postprocessor);
      executioner.run();
    }
  } // end refinement cycle loop
}
