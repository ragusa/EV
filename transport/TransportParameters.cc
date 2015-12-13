/**
 * Constructor.
 */
template<int dim>
TransportParameters<dim>::TransportParameters() :
    problem_id(1),
    end_time(1.0),
    time_step_size(0.001),
    CFL_limit(0.5),
    time_discretization_option(TemporalDiscretization::FE),
    viscosity_option(0),
    entropy_string("0.5*u*u"),
    entropy_derivative_string("u"),
    entropy_residual_coefficient(1.0),
    jump_coefficient(1.0),
    EV_time_discretization(TemporalDiscretization::BE),
    do_not_limit(false),
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

/** \brief defines all of the input parameters
 */
template<int dim>
void TransportParameters<dim>::declare_parameters(ParameterHandler &prm)
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
    prm.declare_entry("Finite element degree", "1", Patterns::Integer(),
        "Polynomial degree of finite elements");
    prm.declare_entry("Number of quadrature points", "3", Patterns::Integer(),
        "Number of quadrature points to use in formula");
  }
  prm.leave_subsection();

  // time parameters
  prm.enter_subsection("time");
  {
    prm.declare_entry("Time discretization option", "FE",
        Patterns::Selection("SS|FE|CN|BE|SSP2|SSP3"), "Choice of temporal discretization");
    prm.declare_entry("End time", "1.0", Patterns::Double(),
        "End time if transient problem is run");
    prm.declare_entry("Time step size", "0.01", Patterns::Double(),
        "Time step size if transient problem is run");
    prm.declare_entry("CFL limit", "0.5", Patterns::Double(),
        "Upper bound for the CFL number");
  }
  prm.leave_subsection();

  // refinement parameters
  prm.enter_subsection("refinement");
  {
    prm.declare_entry("Refinement mode", "space",
        Patterns::Selection("space|time"), "Refinement mode (space or time)");
    prm.declare_entry("Time refinement factor", "0.5", Patterns::Double(),
        "Factor by which to reduce dt if refining time step size");
    prm.declare_entry("Use adaptive refinement", "false", Patterns::Bool(),
        "Option to use adaptive mesh refinement instead of uniform");
    prm.declare_entry("Number of refinement cycles", "5", Patterns::Integer(),
        "Number of mesh refinement cycles");
    prm.declare_entry("Initial refinement level", "2", Patterns::Integer(),
        "Number of refinements for first mesh refinement cycle");
  }
  prm.leave_subsection();

  // linear solver parameters
  prm.enter_subsection("linear solver");
  {
    prm.declare_entry("Linear solver option", "1", Patterns::Integer(),
        "Option for linear solver");
  }
  prm.leave_subsection();

  // nonlinear solver parameters
  prm.enter_subsection("nonlinear solver");
  {
    prm.declare_entry("Nonlinear solver option", "1", Patterns::Integer(),
        "Option for nonlinear solver");
    prm.declare_entry("Nonlinear tolerance", "1.0e-10", Patterns::Double(),
        "Tolerance for nonlinear solver");
    prm.declare_entry("Nonlinear max iteration", "100", Patterns::Integer(),
        "Maximum number of iterations for nonlinear solver");
    prm.declare_entry("Relaxation factor", "1.0", Patterns::Double(),
        "Relaxation factor for nonlinear solver updates");
  }
  prm.leave_subsection();

  // viscosity parameters
  prm.enter_subsection("viscosity");
  {
    prm.declare_entry("Viscosity option", "0", Patterns::Integer(),
        "Option for viscosity to be used");
  }
  prm.leave_subsection();
  // viscosity parameters
  prm.enter_subsection("viscosity");
  {
    prm.declare_entry("Viscosity option", "0", Patterns::Integer(),
        "Option for viscosity to be used");
  }
  prm.leave_subsection();

  // entropy viscosity parameters
  prm.enter_subsection("entropy viscosity");
  {
    prm.declare_entry("Entropy string", "0.5*u*u", Patterns::Anything(),
        "String for entropy function");
    prm.declare_entry("Entropy derivative string", "0.5*u*u",
        Patterns::Anything(), "String for entropy derivative function");
    prm.declare_entry("Entropy residual coefficient", "1.0", Patterns::Double(),
        "Coefficient for the entropy viscosity");
    prm.declare_entry("Jump coefficient", "1.0", Patterns::Double(),
        "Coefficient for jumps used with entropy viscosity");
    prm.declare_entry("EV temporal discretization", "BE",
        Patterns::Selection("FE|BE|CN|BDF2"),
        "Temporal discretization for entropy residual");
  }
  prm.leave_subsection();

  // output parameters
  prm.enter_subsection("output");
  {
    prm.declare_entry("Output mesh", "false", Patterns::Bool(),
        "Option to output meshes as .eps files");
    prm.declare_entry("Output exact solution", "false", Patterns::Bool(),
        "Option to output exact solution");
    prm.declare_entry("Exact solution refinement level", "5", Patterns::Integer(),
        "Refinement level for exact solution");
    prm.declare_entry("Output initial solution", "false", Patterns::Bool(),
        "Option to output initial solution");
    prm.declare_entry("Output DMP bounds", "false", Patterns::Bool(),
        "Option to output DMP bounds");
    prm.declare_entry("Save convergence results", "false", Patterns::Bool(),
        "Option to save convegence results");
    prm.declare_entry("Print solution", "false", Patterns::Bool(),
        "Option to print final solution");
  }
  prm.leave_subsection();

  // fct parameters
  prm.enter_subsection("fct");
  {
    prm.declare_entry("Do not limit", "false", Patterns::Bool(),
        "Option to choose not to limit when constructing high-order solution");
  }
  prm.leave_subsection();
}

/** \brief get the input parameters
 */
template<int dim>
void TransportParameters<dim>::get_parameters(ParameterHandler &prm)
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
    std::string time_discretization_string = prm.get("Time discretization option");
    if (time_discretization_string == "SS")
    {
      time_discretization_option = TemporalDiscretization::SS;
    }
    else if (time_discretization_string == "FE")
    {
      time_discretization_option = TemporalDiscretization::FE;
    }
    else if (time_discretization_string == "CN")
    {
      time_discretization_option = TemporalDiscretization::CN;
    }
    else if (time_discretization_string == "BE")
    {
      time_discretization_option = TemporalDiscretization::BE;
    }
    else if (time_discretization_string == "SSP2")
    {
      time_discretization_option = TemporalDiscretization::SSP2;
    }
    else if (time_discretization_string == "SSP3")
    {
      time_discretization_option = TemporalDiscretization::SSP3;
    }
    else
    {
      ExcNotImplemented();
    }
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
      ExcNotImplemented();
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

    std::string EV_time_discretization_string = prm.get(
        "EV temporal discretization");
    if (EV_time_discretization_string == "FE")
    {
      EV_time_discretization = TemporalDiscretization::FE;
    }
    else if (EV_time_discretization_string == "BE")
    {
      EV_time_discretization = TemporalDiscretization::BE;
    }
    else if (EV_time_discretization_string == "CN")
    {
      EV_time_discretization = TemporalDiscretization::CN;
    }
    else if (EV_time_discretization_string == "BDF2")
    {
      EV_time_discretization = TemporalDiscretization::BDF2;
    }
    else
    {
      ExcNotImplemented();
    }
  }
  prm.leave_subsection();

  // output parameters
  prm.enter_subsection("output");
  {
    output_mesh = prm.get_bool("Output mesh");
    output_exact_solution = prm.get_bool("Output exact solution");
    exact_solution_refinement_level = prm.get_integer(
        "Exact solution refinement level");
    output_initial_solution = prm.get_bool("Output initial solution");
    output_DMP_bounds = prm.get_bool("Output DMP bounds");
    save_convergence_results = prm.get_bool("Save convergence results");
    print_solution = prm.get_bool("Print solution");
  }
  prm.leave_subsection();

  // fct parameters
  prm.enter_subsection("fct");
  {
    do_not_limit = prm.get_bool("Do not limit");
  }
  prm.leave_subsection();
}
