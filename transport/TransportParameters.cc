/**
 * Constructor.
 */
template <int dim>
TransportParameters<dim>::TransportParameters()
  : temporal_discretization(TemporalDiscretization::ssprk),
    entropy_temporal_discretization(EntropyTemporalDiscretization::BE),
    ssprk_method(SSPRKMethod::FE),
    theta_method(ThetaMethod::FE),
    theta(0.0),
    fct_initialization_option(FCTInitializationOption::zero),
    problem_id(1),
    end_time(1.0),
    time_step_size(0.001),
    CFL_limit(0.5),
    viscosity_option(0),
    entropy_string("0.5*u*u"),
    entropy_derivative_string("u"),
    entropy_residual_coefficient(1.0),
    jump_coefficient(1.0),
    do_not_limit(false),
    skip_fct_if_bounds_satisfied(false),
    use_cumulative_antidiffusion_algorithm(false),
    include_analytic_bounds(false),
    refinement_mode(RefinementMode::space),
    time_refinement_factor(0.5),
    use_adaptive_refinement(false),
    initial_refinement_level(2),
    n_refinement_cycles(5),
    degree(1),
    n_quadrature_points(3),
    linear_solver_option(1),
    nonlinear_solver_option(1),
    nonlinear_tolerance(1.0e-10),
    nonlinear_max_iteration(100),
    relaxation_factor(1.0),
    output_mesh(false),
    output_exact_solution(false),
    exact_solution_refinement_level(5),
    output_initial_solution(false),
    output_DMP_bounds(false),
    save_convergence_results(false),
    print_solution(false)
{
}

/**
 * \brief defines all of the input parameters
 */
template <int dim>
void TransportParameters<dim>::declare_parameters(ParameterHandler & prm)
{
  // problem parameters
  prm.enter_subsection("problem");
  {
    prm.declare_entry("Problem ID", "1", Patterns::Integer(), "Problem ID");
  }
  prm.leave_subsection();

  // finite element parameters
  prm.enter_subsection("finite element");
  {
    prm.declare_entry("Finite element degree",
                      "1",
                      Patterns::Integer(),
                      "Polynomial degree of finite elements");
    prm.declare_entry("Number of quadrature points",
                      "3",
                      Patterns::Integer(),
                      "Number of quadrature points to use in formula");
  }
  prm.leave_subsection();

  // time parameters
  prm.enter_subsection("time");
  {
    prm.declare_entry("Time discretization",
                      "ssprk",
                      Patterns::Selection("ss|theta|ssprk"),
                      "Choice of temporal discretization");
    prm.declare_entry("Entropy time discretization",
                      "BE",
                      Patterns::Selection("FE|CN|BE|BDF2"),
                      "Choice of temporal discretization for entropy");
    prm.declare_entry("SSPRK method",
                      "FE",
                      Patterns::Selection("FE|SSP2|SSP3"),
                      "Choice of SSPRK method");
    prm.declare_entry("Theta method",
                      "FE",
                      Patterns::Selection("FE|CN|BE"),
                      "Choice of theta time discretization method");
    prm.declare_entry("End time",
                      "1.0",
                      Patterns::Double(),
                      "End time if transient problem is run");
    prm.declare_entry("Time step size",
                      "0.01",
                      Patterns::Double(),
                      "Time step size if transient problem is run");
    prm.declare_entry(
      "CFL limit", "0.5", Patterns::Double(), "Upper bound for the CFL number");
  }
  prm.leave_subsection();

  // refinement parameters
  prm.enter_subsection("refinement");
  {
    prm.declare_entry("Refinement mode",
                      "space",
                      Patterns::Selection("space|time"),
                      "Refinement mode (space or time)");
    prm.declare_entry("Time refinement factor",
                      "0.5",
                      Patterns::Double(),
                      "Factor by which to reduce dt if refining time step size");
    prm.declare_entry(
      "Use adaptive refinement",
      "false",
      Patterns::Bool(),
      "Option to use adaptive mesh refinement instead of uniform");
    prm.declare_entry("Number of refinement cycles",
                      "5",
                      Patterns::Integer(),
                      "Number of mesh refinement cycles");
    prm.declare_entry("Initial refinement level",
                      "2",
                      Patterns::Integer(),
                      "Number of refinements for first mesh refinement cycle");
  }
  prm.leave_subsection();

  // linear solver parameters
  prm.enter_subsection("linear solver");
  {
    prm.declare_entry("Linear solver option",
                      "1",
                      Patterns::Integer(),
                      "Option for linear solver");
  }
  prm.leave_subsection();

  // nonlinear solver parameters
  prm.enter_subsection("nonlinear solver");
  {
    prm.declare_entry("Nonlinear solver option",
                      "1",
                      Patterns::Integer(),
                      "Option for nonlinear solver");
    prm.declare_entry("Nonlinear tolerance",
                      "1.0e-10",
                      Patterns::Double(),
                      "Tolerance for nonlinear solver");
    prm.declare_entry("Nonlinear max iteration",
                      "100",
                      Patterns::Integer(),
                      "Maximum number of iterations for nonlinear solver");
    prm.declare_entry("Relaxation factor",
                      "1.0",
                      Patterns::Double(),
                      "Relaxation factor for nonlinear solver updates");
  }
  prm.leave_subsection();

  // viscosity parameters
  prm.enter_subsection("viscosity");
  {
    prm.declare_entry("Viscosity option",
                      "0",
                      Patterns::Integer(),
                      "Option for viscosity to be used");
  }
  prm.leave_subsection();
  // viscosity parameters
  prm.enter_subsection("viscosity");
  {
    prm.declare_entry("Viscosity option",
                      "0",
                      Patterns::Integer(),
                      "Option for viscosity to be used");
  }
  prm.leave_subsection();

  // entropy viscosity parameters
  prm.enter_subsection("entropy viscosity");
  {
    prm.declare_entry("Entropy string",
                      "0.5*u*u",
                      Patterns::Anything(),
                      "String for entropy function");
    prm.declare_entry("Entropy derivative string",
                      "0.5*u*u",
                      Patterns::Anything(),
                      "String for entropy derivative function");
    prm.declare_entry("Entropy residual coefficient",
                      "1.0",
                      Patterns::Double(),
                      "Coefficient for the entropy viscosity");
    prm.declare_entry("Jump coefficient",
                      "1.0",
                      Patterns::Double(),
                      "Coefficient for jumps used with entropy viscosity");
    prm.declare_entry("EV temporal discretization",
                      "BE",
                      Patterns::Selection("BE|CN|BDF2"),
                      "Temporal discretization for entropy residual");
  }
  prm.leave_subsection();

  // output parameters
  prm.enter_subsection("output");
  {
    prm.declare_entry("Output mesh",
                      "false",
                      Patterns::Bool(),
                      "Option to output meshes as .eps files");
    prm.declare_entry("Output exact solution",
                      "false",
                      Patterns::Bool(),
                      "Option to output exact solution");
    prm.declare_entry("Exact solution refinement level",
                      "5",
                      Patterns::Integer(),
                      "Refinement level for exact solution");
    prm.declare_entry("Output initial solution",
                      "false",
                      Patterns::Bool(),
                      "Option to output initial solution");
    prm.declare_entry("Output DMP bounds",
                      "false",
                      Patterns::Bool(),
                      "Option to output DMP bounds");
    prm.declare_entry("Save convergence results",
                      "false",
                      Patterns::Bool(),
                      "Option to save convegence results");
    prm.declare_entry("Print solution",
                      "false",
                      Patterns::Bool(),
                      "Option to print final solution");
  }
  prm.leave_subsection();

  // fct parameters
  prm.enter_subsection("fct");
  {
    prm.declare_entry("FCT initialization option",
                      "zero",
                      Patterns::Anything(),
                      "Initialization option for implicit FCT schemes");
    prm.declare_entry(
      "Do not limit",
      "false",
      Patterns::Bool(),
      "Option to choose not to limit when constructing high-order solution");
    prm.declare_entry(
      "Skip fct if bounds satisfied",
      "false",
      Patterns::Bool(),
      "Option to skip FCT if high-order solution satisfies bounds");
    prm.declare_entry(
      "Use cumulative antidiffusion algorithm",
      "false",
      Patterns::Bool(),
      "Option to use cumulative antidiffusion algorithm for implicit FCT");
    prm.declare_entry("Include analytic bounds",
                      "false",
                      Patterns::Bool(),
                      "Option to extend FCT bounds using analytic DMP");
  }
  prm.leave_subsection();
}

/** \brief get the input parameters
 */
template <int dim>
void TransportParameters<dim>::get_parameters(ParameterHandler & prm)
{
  // problem parameters
  prm.enter_subsection("problem");
  {
    problem_id = prm.get_integer("Problem ID");
  }
  prm.leave_subsection();

  // finite element parameters
  prm.enter_subsection("finite element");
  {
    degree = prm.get_integer("Finite element degree");
    n_quadrature_points = prm.get_integer("Number of quadrature points");
  }
  prm.leave_subsection();

  // time parameters
  prm.enter_subsection("time");
  {
    // get temporal discretization
    std::string temporal_discretization_string = prm.get("Time discretization");
    if (temporal_discretization_string == "ss")
    {
      temporal_discretization = TemporalDiscretization::ss;
    }
    else if (temporal_discretization_string == "theta")
    {
      temporal_discretization = TemporalDiscretization::theta;
    }
    else if (temporal_discretization_string == "ssprk")
    {
      temporal_discretization = TemporalDiscretization::ssprk;
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }

    // get temporal discretization for entropy
    std::string entropy_temporal_discretization_string =
      prm.get("Entropy time discretization");
    if (entropy_temporal_discretization_string == "CN")
    {
      entropy_temporal_discretization = EntropyTemporalDiscretization::CN;
    }
    else if (entropy_temporal_discretization_string == "BE")
    {
      entropy_temporal_discretization = EntropyTemporalDiscretization::BE;
    }
    else if (entropy_temporal_discretization_string == "BDF2")
    {
      entropy_temporal_discretization = EntropyTemporalDiscretization::BDF2;
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }

    // get SSPRK method
    std::string ssprk_method_string = prm.get("SSPRK method");
    if (ssprk_method_string == "FE")
    {
      ssprk_method = SSPRKMethod::FE;
    }
    else if (ssprk_method_string == "SSP2")
    {
      ssprk_method = SSPRKMethod::SSP2;
    }
    else if (ssprk_method_string == "SSP3")
    {
      ssprk_method = SSPRKMethod::SSP3;
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }

    // get theta method and corresponding parameter
    std::string theta_method_string = prm.get("Theta method");
    if (theta_method_string == "FE")
    {
      theta_method = ThetaMethod::FE;
      theta = 0.0;
    }
    else if (theta_method_string == "CN")
    {
      theta_method = ThetaMethod::CN;
      theta = 0.5;
    }
    else if (theta_method_string == "BE")
    {
      theta_method = ThetaMethod::BE;
      theta = 1.0;
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }

    // get other parameters
    end_time = prm.get_double("End time");
    time_step_size = prm.get_double("Time step size");
    CFL_limit = prm.get_double("CFL limit");
  }
  prm.leave_subsection();

  // refinement parameters
  prm.enter_subsection("refinement");
  {
    std::string refinement_mode_string = prm.get("Refinement mode");
    if (refinement_mode_string == "space")
    {
      refinement_mode = RefinementMode::space;
    }
    else if (refinement_mode_string == "time")
    {
      refinement_mode = RefinementMode::time;
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }

    time_refinement_factor = prm.get_double("Time refinement factor");
    use_adaptive_refinement = prm.get_bool("Use adaptive refinement");
    n_refinement_cycles = prm.get_integer("Number of refinement cycles");
    initial_refinement_level = prm.get_integer("Initial refinement level");
  }
  prm.leave_subsection();

  // linear solver parameters
  prm.enter_subsection("linear solver");
  {
    linear_solver_option = prm.get_integer("Linear solver option");
  }
  prm.leave_subsection();

  // nonlinear solver parameters
  prm.enter_subsection("nonlinear solver");
  {
    nonlinear_solver_option = prm.get_integer("Nonlinear solver option");
    nonlinear_tolerance = prm.get_double("Nonlinear tolerance");
    nonlinear_max_iteration = prm.get_integer("Nonlinear max iteration");
    relaxation_factor = prm.get_double("Relaxation factor");
  }
  prm.leave_subsection();

  // viscosity parameters
  prm.enter_subsection("viscosity");
  {
    viscosity_option = prm.get_integer("Viscosity option");
  }
  prm.leave_subsection();

  // entropy viscosity parameters
  prm.enter_subsection("entropy viscosity");
  {
    entropy_string = prm.get("Entropy string");
    entropy_derivative_string = prm.get("Entropy derivative string");
    entropy_residual_coefficient = prm.get_double("Entropy residual coefficient");
    jump_coefficient = prm.get_double("Jump coefficient");
  }
  prm.leave_subsection();

  // output parameters
  prm.enter_subsection("output");
  {
    output_mesh = prm.get_bool("Output mesh");
    output_exact_solution = prm.get_bool("Output exact solution");
    exact_solution_refinement_level =
      prm.get_integer("Exact solution refinement level");
    output_initial_solution = prm.get_bool("Output initial solution");
    output_DMP_bounds = prm.get_bool("Output DMP bounds");
    save_convergence_results = prm.get_bool("Save convergence results");
    print_solution = prm.get_bool("Print solution");
  }
  prm.leave_subsection();

  // fct parameters
  prm.enter_subsection("fct");
  {
    std::string fct_initialization_string = prm.get("FCT initialization option");
    if (fct_initialization_string == "zero")
    {
      fct_initialization_option = FCTInitializationOption::zero;
    }
    else if (fct_initialization_string == "low")
    {
      fct_initialization_option = FCTInitializationOption::low;
    }
    else if (fct_initialization_string == "high")
    {
      fct_initialization_option = FCTInitializationOption::high;
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }

    do_not_limit = prm.get_bool("Do not limit");

    skip_fct_if_bounds_satisfied = prm.get_bool("Skip fct if bounds satisfied");

    use_cumulative_antidiffusion_algorithm =
      prm.get_bool("Use cumulative antidiffusion algorithm");

    include_analytic_bounds = prm.get_bool("Include analytic bounds");
  }
  prm.leave_subsection();
}
