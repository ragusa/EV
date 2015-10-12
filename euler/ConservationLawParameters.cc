/**
 * \file ConservationLawParameters.cc
 * \brief Provides the function definitions for the ConservationLawParameters
 *        class.
 */

/**
 *  \brief Declares the parameters that will be in input file.
 *
 *  This function declares all of the input parameters required
 *  for the ConservationLaw class.
 *  \param prm parameter handler for conservation law parameters
 */
template <int dim>
void ConservationLawParameters<dim>::declare_conservation_law_parameters(
  ParameterHandler & prm)
{
  // finite element parameters
  prm.enter_subsection("finite element");
  {
    prm.declare_entry(
      "degree", "1", Patterns::Integer(), "polynomial degree of finite elements");
  }
  prm.leave_subsection();

  // quadrature
  prm.enter_subsection("quadrature");
  {
    prm.declare_entry("number of quadrature points",
                      "3",
                      Patterns::Integer(),
                      "number of quadrature points per dimension");
  }
  prm.leave_subsection();

  // refinement parameters
  prm.enter_subsection("refinement");
  {
    prm.declare_entry("Refinement mode",
                      "space",
                      Patterns::Selection("space|time"),
                      "Refinement mode (space or time)");

    prm.declare_entry("initial refinement level",
                      "3",
                      Patterns::Integer(),
                      "initial number of uniform refinements");
    prm.declare_entry(
      "refinement cycles",
      "2",
      Patterns::Integer(),
      "Number of refinement cycles. Must be >= 1. '1' Corresponds to"
      " having only the initial uniform refinement cycles.");
    prm.declare_entry(
      "refinement fraction",
      "0.3",
      Patterns::Double(),
      "Fraction of cells to refine in an adaptive refinement step");
    prm.declare_entry(
      "coarsening fraction",
      "0.03",
      Patterns::Double(),
      "Fraction of cells to coarsen in an adaptive refinement step");
  }
  prm.leave_subsection();

  // time parameters
  prm.enter_subsection("time");
  {
    prm.declare_entry("time step size method",
                      "cfl_condition",
                      Patterns::Selection("constant_dt|cfl_condition"),
                      "method of computing time step size");
    prm.declare_entry("use default end time",
                      "true",
                      Patterns::Bool(),
                      "Use default end time for the particular test problem");
    prm.declare_entry(
      "final time", "1.0", Patterns::Double(), "final time value");
    prm.declare_entry(
      "time step size", "1e-3", Patterns::Double(), "time step size");
    prm.declare_entry(
      "cfl", "0.5", Patterns::Double(), "CFL number to be used if CFL condition "
                                        "is used to compute time step size.");
    prm.declare_entry("steady state tolerance",
                      "1.0e-6",
                      Patterns::Double(),
                      "Tolerance for achievement of steady-state");
  }
  prm.leave_subsection();

  // temporal integrator
  prm.enter_subsection("temporal integration");
  {
    prm.declare_entry("temporal integrator",
                      "runge_kutta",
                      Patterns::Selection("runge_kutta"),
                      "The method used for advancing time. "
                      "Choices are <runge_kutta>.");
    prm.declare_entry("runge kutta method",
                      "ERK1",
                      Patterns::Selection("ERK1|ERK2|ERK3|ERK4|SDIRK22"),
                      "Runge-Kutta method to use.");
  }
  prm.leave_subsection();

  // nonlinear solver parameters
  prm.enter_subsection("nonlinear solver");
  {
    prm.declare_entry(
      "nonlinear verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose"),
      "State whether output from nonlinear solver runs should be printed. "
      "Choices are <quiet|verbose>.");
    prm.declare_entry("nonlinear method",
                      "newton",
                      Patterns::Selection("newton"),
                      "The kind of nonlinear solver for the linear system. "
                      "Choices are <newton>.");
    prm.declare_entry("nonlinear absolute tolerance",
                      "1e-10",
                      Patterns::Double(),
                      "Nonlinear absolute tolerance");
    prm.declare_entry("nonlinear relative tolerance",
                      "1e-10",
                      Patterns::Double(),
                      "Nonlinear relative tolerance");
    prm.declare_entry("max nonlinear iterations",
                      "300",
                      Patterns::Integer(),
                      "Maximum nonlinear iterations");
    prm.declare_entry("damping", "1.0", Patterns::Double(), "damping");
  }
  prm.leave_subsection();

  // linear solver parameters
  prm.enter_subsection("linear solver");
  {
    prm.declare_entry(
      "linear verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose"),
      "State whether output from linear solver runs should be printed. "
      "Choices are <quiet|verbose>.");
    prm.declare_entry("linear method",
                      "direct",
                      Patterns::Selection("direct|gmres|cg"),
                      "The kind of linear solver for the linear system. "
                      "Choices are <direct|gmres|cg>.");
    prm.declare_entry(
      "mass matrix linear method",
      "cg",
      Patterns::Selection("direct|gmres|cg"),
      "The linear solver used to implicitly invert the mass matrix. "
      "Choices are <direct|gmres|cg>.");
    prm.declare_entry("linear absolute tolerance",
                      "1e-10",
                      Patterns::Double(),
                      "Linear absolute tolerance");
    prm.declare_entry("linear relative tolerance",
                      "1e-10",
                      Patterns::Double(),
                      "Linear relative tolerance");
    prm.declare_entry("max linear iterations",
                      "300",
                      Patterns::Integer(),
                      "Maximum linear solver iterations");
  }
  prm.leave_subsection();

  // artificial viscosity
  prm.enter_subsection("artificial viscosity");
  {
    prm.declare_entry("viscosity type",
                      "constant",
                      Patterns::Anything(),
                      "choice for artificial viscosity");
    prm.declare_entry("constant viscosity value",
                      "1e-3",
                      Patterns::Double(),
                      "viscosity value if constant viscosity chosen");
    prm.declare_entry(
      "first order viscosity coefficient",
      "1e-3",
      Patterns::Double(),
      "tuning constant value to be used with first-order viscosity");
    prm.declare_entry("entropy viscosity coefficient",
                      "1e-3",
                      Patterns::Double(),
                      "tuning constant value to be used with entropy viscosity");
    prm.declare_entry("jump coefficient",
                      "1.0",
                      Patterns::Double(),
                      "tuning constant value to be used with jumps");
  }
  prm.leave_subsection();

  // output
  prm.enter_subsection("output");
  {
    prm.declare_entry("output period",
                      "1",
                      Patterns::Integer(),
                      "Period of time steps for outputting the solution, e.g.,"
                      " 1 would output every time step,"
                      " and 2 would output every other time step, etc.");
    prm.declare_entry("output mesh",
                      "true",
                      Patterns::Bool(),
                      "option to output the mesh to a file");
    prm.declare_entry("output mass matrix",
                      "true",
                      Patterns::Bool(),
                      "option to output the mass matrix to a file");
    prm.declare_entry("output viscosity",
                      "true",
                      Patterns::Bool(),
                      "option to output the viscosity to a file");
    prm.declare_entry("output exact solution",
                      "true",
                      Patterns::Bool(),
                      "option to output the exact solution if it exists");
    prm.declare_entry("save convergence results",
                      "true",
                      Patterns::Bool(),
                      "option to save convergence results to a file");
    prm.declare_entry("exact solution refinement level",
                      "7",
                      Patterns::Integer(),
                      "refinement level to be used for the exact solution");
  }
  prm.leave_subsection();
}

/** \fn
 * ConservationLawParameters<dim>::get_conservation_law_parameters(ParameterHandler
 * &prm)
 *  \brief Gets the parameters from the parameter handler.
 *
 *  This function takes the input parameters from the parameter
 *  handler into the member variables.
 *  \param prm parameter handler for conservation law parameters
 */
template <int dim>
void ConservationLawParameters<dim>::get_conservation_law_parameters(
  ParameterHandler & prm)
{
  // finite element parameters
  prm.enter_subsection("finite element");
  {
    degree = prm.get_integer("degree");
  }
  prm.leave_subsection();

  // quadrature
  prm.enter_subsection("quadrature");
  {
    n_quadrature_points = prm.get_integer("number of quadrature points");
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

    initial_refinement_level = prm.get_integer("initial refinement level");
    n_refinement_cycles = prm.get_integer("refinement cycles");
    refinement_fraction = prm.get_double("refinement fraction");
    coarsening_fraction = prm.get_double("coarsening fraction");
  }
  prm.leave_subsection();

  // time parameters
  prm.enter_subsection("time");
  {
    const std::string time_choice = prm.get("time step size method");
    if (time_choice == "constant_dt")
      time_step_size_method = constant_dt;
    else if (time_choice == "cfl_condition")
      time_step_size_method = cfl_condition;
    else
      Assert(false, ExcNotImplemented());

    use_default_end_time = prm.get_bool("use default end time");
    end_time = prm.get_double("final time");
    time_step_size = prm.get_double("time step size");
    cfl = prm.get_double("cfl");
    steady_state_tolerance = prm.get_double("steady state tolerance");
  }
  prm.leave_subsection();

  // temporal integrator
  prm.enter_subsection("temporal integration");
  {
    const std::string temporal_choice = prm.get("temporal integrator");
    if (temporal_choice == "runge_kutta")
      temporal_integrator = runge_kutta;
    else
      Assert(false, ExcNotImplemented());

    const std::string rk_choice = prm.get("runge kutta method");
    if (rk_choice == "ERK1")
      time_discretization = ERK1;
    else if (rk_choice == "ERK2")
      time_discretization = ERK2;
    else if (rk_choice == "ERK3")
      time_discretization = ERK3;
    else if (rk_choice == "ERK4")
      time_discretization = ERK4;
    else if (rk_choice == "SDIRK22")
      time_discretization = SDIRK22;
    else
      Assert(false, ExcNotImplemented());
  }
  prm.leave_subsection();

  // nonlinear solver parameters
  prm.enter_subsection("nonlinear solver");
  {
    const std::string verbosity = prm.get("nonlinear verbosity");
    if (verbosity == "verbose")
      nonlinear_verbosity = verbose;
    if (verbosity == "quiet")
      nonlinear_verbosity = quiet;

    const std::string solver = prm.get("nonlinear method");
    if (solver == "newton")
      nonlinear_solver = newton;

    nonlinear_atol = prm.get_double("nonlinear absolute tolerance");
    nonlinear_rtol = prm.get_double("nonlinear relative tolerance");
    max_nonlinear_iterations = prm.get_integer("max nonlinear iterations");
    damping = prm.get_double("damping");
  }
  prm.leave_subsection();

  // linear solver parameters
  prm.enter_subsection("linear solver");
  {
    const std::string verbosity = prm.get("linear verbosity");
    if (verbosity == "verbose")
      linear_verbosity = verbose;
    if (verbosity == "quiet")
      linear_verbosity = quiet;

    const std::string solver = prm.get("linear method");
    if (solver == "direct")
      linear_solver = direct;
    else if (solver == "gmres")
      linear_solver = gmres;
    else if (solver == "cg")
      linear_solver = cg;

    const std::string mass_solver = prm.get("mass matrix linear method");
    if (mass_solver == "direct")
      mass_matrix_linear_solver = direct;
    else if (mass_solver == "gmres")
      mass_matrix_linear_solver = gmres;
    else if (mass_solver == "cg")
      mass_matrix_linear_solver = cg;

    linear_atol = prm.get_double("linear absolute tolerance");
    linear_rtol = prm.get_double("linear relative tolerance");
    max_linear_iterations = prm.get_integer("max linear iterations");
  }
  prm.leave_subsection();

  // artificial viscosity
  prm.enter_subsection("artificial viscosity");
  {
    const std::string viscosity_choice = prm.get("viscosity type");
    if (viscosity_choice == "none")
      viscosity_type = none;
    else if (viscosity_choice == "constant")
      viscosity_type = constant;
    else if (viscosity_choice == "old_first_order")
      viscosity_type = old_first_order;
    else if (viscosity_choice == "max_principle")
      viscosity_type = max_principle;
    else if (viscosity_choice == "entropy")
      viscosity_type = entropy;
    else
      Assert(false, ExcNotImplemented());

    constant_viscosity_value = prm.get_double("constant viscosity value");
    first_order_viscosity_coef =
      prm.get_double("first order viscosity coefficient");
    entropy_viscosity_coef = prm.get_double("entropy viscosity coefficient");
    jump_coef = prm.get_double("jump coefficient");
  }
  prm.leave_subsection();

  // output parameters
  prm.enter_subsection("output");
  {
    output_period = prm.get_integer("output period");
    output_mesh = prm.get_bool("output mesh");
    output_mass_matrix = prm.get_bool("output mass matrix");
    output_viscosity = prm.get_bool("output viscosity");
    output_exact_solution = prm.get_bool("output exact solution");
    exact_solution_refinement_level =
      prm.get_integer("exact solution refinement level");
    save_convergence_results = prm.get_bool("save convergence results");
  }
  prm.leave_subsection();
}
