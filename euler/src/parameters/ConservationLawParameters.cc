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
  // scheme
  prm.enter_subsection("scheme");
  {
    prm.declare_entry(
      "scheme", "high", Patterns::Selection("low|high|fct"), "scheme");
    prm.declare_entry(
      "low order scheme",
      "standard",
      Patterns::Selection("constant|standard|dmp|di_visc|di_diff"),
      "low-order scheme");
    prm.declare_entry("high order scheme",
                      "galerkin",
                      Patterns::Selection("galerkin|entropy_visc|entropy_diff"),
                      "high-order scheme");
  }
  prm.leave_subsection();

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
                      "constant",
                      Patterns::Selection("constant|cfl"),
                      "method of computing time step size");
    prm.declare_entry("use default end time",
                      "true",
                      Patterns::Bool(),
                      "Use default end time for the particular test problem");
    prm.declare_entry(
      "final time", "1.0", Patterns::Double(), "final time value");
    prm.declare_entry(
      "time step size", "1e-3", Patterns::Double(), "time step size");
    prm.declare_entry("cfl",
                      "0.5",
                      Patterns::Double(),
                      "CFL number to be used if CFL condition "
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
                      "FE",
                      Patterns::Selection("FE|SSP2|SSP3"),
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
    prm.declare_entry("constant viscosity value",
                      "1e-3",
                      Patterns::Double(),
                      "viscosity value if constant viscosity chosen");
    prm.declare_entry(
      "first order viscosity coefficient",
      "1e-3",
      Patterns::Double(),
      "tuning constant value to be used with first-order viscosity");
    prm.declare_entry("use low order viscosity for first time step",
                      "true",
                      Patterns::Bool(),
                      "option to use low-order viscosity for first time step"
                      " of entropy viscosity method");
    prm.declare_entry("entropy residual coefficient",
                      "1.0",
                      Patterns::Double(),
                      "tuning constant value to be used with entropy residual");
    prm.declare_entry("entropy jump coefficient",
                      "1.0",
                      Patterns::Double(),
                      "tuning constant value to be used with entropy jumps");
    prm.declare_entry("entropy viscosity smoothing",
                      "none",
                      Patterns::Anything(),
                      "Type of smoothing to apply to entropy viscosity");
    prm.declare_entry(
      "entropy viscosity smoothing weight",
      "2",
      Patterns::Integer(),
      "Weight for center of Laplacian smoothing for entropy viscosity");
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
    prm.declare_entry("max transient output size",
                      "1e9",
                      Patterns::Double(),
                      "Maximum number of bytes to use for transient output"
                      " files");
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
    prm.declare_entry("output viscosity transient",
                      "false",
                      Patterns::Bool(),
                      "option to output the viscosity transient");
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
    prm.declare_entry("print final solution",
                      "false",
                      Patterns::Bool(),
                      "option to print final solution");
    prm.declare_entry(
      "verbosity level", "1", Patterns::Integer(), "level of verbosity");
  }
  prm.leave_subsection();

  // fct
  prm.enter_subsection("fct");
  {
    prm.declare_entry("antidiffusion",
                      "limited",
                      Patterns::Selection("limited|full|none"),
                      "Option for antidiffusion in FCT scheme");
    prm.declare_entry(
      "fct synchronization",
      "none",
      Patterns::Selection("none|min|compound"),
      "Option for synchronization of limiting coefficients in FCT scheme");
    prm.declare_entry("use star states in fct bounds",
                      "false",
                      Patterns::Bool(),
                      "Option to include star states in FCT bounds");
  }
  prm.leave_subsection();
}

/**
 * \brief Gets the parameters from the parameter handler.
 *
 * This function takes the input parameters from the parameter
 * handler into the member variables.
 *
 * \param[in] prm parameter handler for conservation law parameters
 */
template <int dim>
void ConservationLawParameters<dim>::get_conservation_law_parameters(
  ParameterHandler & prm)
{
  // scheme
  prm.enter_subsection("scheme");
  {
    // scheme
    std::string scheme_string = prm.get("scheme");
    if (scheme_string == "low")
      scheme = Scheme::low;
    else if (scheme_string == "high")
      scheme = Scheme::high;
    else if (scheme_string == "fct")
      scheme = Scheme::fct;
    else
      Assert(false, ExcNotImplemented());

    // low-order scheme
    std::string low_order_scheme_string = prm.get("low order scheme");
    if (low_order_scheme_string == "constant")
      low_order_scheme = LowOrderScheme::constant;
    else if (low_order_scheme_string == "standard")
      low_order_scheme = LowOrderScheme::standard;
    else if (low_order_scheme_string == "dmp")
      low_order_scheme = LowOrderScheme::dmp;
    else if (low_order_scheme_string == "di_visc")
      low_order_scheme = LowOrderScheme::di_visc;
    else if (low_order_scheme_string == "di_diff")
      low_order_scheme = LowOrderScheme::di_diff;
    else
      Assert(false, ExcNotImplemented());

    // high-order scheme
    std::string high_order_scheme_string = prm.get("high order scheme");
    if (high_order_scheme_string == "galerkin")
      high_order_scheme = HighOrderScheme::galerkin;
    else if (high_order_scheme_string == "entropy_visc")
      high_order_scheme = HighOrderScheme::entropy_visc;
    else if (high_order_scheme_string == "entropy_diff")
      high_order_scheme = HighOrderScheme::entropy_diff;
    else
      Assert(false, ExcNotImplemented());
  }
  prm.leave_subsection();

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
    if (time_choice == "constant")
      time_step_size_method = TimeStepSizeMethod::constant;
    else if (time_choice == "cfl")
      time_step_size_method = TimeStepSizeMethod::cfl;
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
      temporal_integrator = TemporalIntegrator::runge_kutta;
    else
      Assert(false, ExcNotImplemented());

    const std::string rk_choice = prm.get("runge kutta method");
    if (rk_choice == "FE")
      time_discretization = TemporalDiscretization::FE;
    else if (rk_choice == "SSP2")
      time_discretization = TemporalDiscretization::SSP2;
    else if (rk_choice == "SSP3")
      time_discretization = TemporalDiscretization::SSP3;
    else
      Assert(false, ExcNotImplemented());
  }
  prm.leave_subsection();

  // nonlinear solver parameters
  prm.enter_subsection("nonlinear solver");
  {
    const std::string verbosity = prm.get("nonlinear verbosity");
    if (verbosity == "verbose")
      nonlinear_verbosity = Verbosity::verbose;
    if (verbosity == "quiet")
      nonlinear_verbosity = Verbosity::quiet;

    const std::string solver = prm.get("nonlinear method");
    if (solver == "newton")
      nonlinear_solver = NonlinearSolverType::newton;

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
      linear_verbosity = Verbosity::verbose;
    if (verbosity == "quiet")
      linear_verbosity = Verbosity::quiet;

    const std::string solver = prm.get("linear method");
    if (solver == "direct")
      linear_solver = LinearSolverType::direct;
    else if (solver == "gmres")
      linear_solver = LinearSolverType::gmres;
    else if (solver == "cg")
      linear_solver = LinearSolverType::cg;

    const std::string mass_solver = prm.get("mass matrix linear method");
    if (mass_solver == "direct")
      mass_matrix_linear_solver = LinearSolverType::direct;
    else if (mass_solver == "gmres")
      mass_matrix_linear_solver = LinearSolverType::gmres;
    else if (mass_solver == "cg")
      mass_matrix_linear_solver = LinearSolverType::cg;

    linear_atol = prm.get_double("linear absolute tolerance");
    linear_rtol = prm.get_double("linear relative tolerance");
    max_linear_iterations = prm.get_integer("max linear iterations");
  }
  prm.leave_subsection();

  // artificial viscosity
  prm.enter_subsection("artificial viscosity");
  {
    constant_viscosity_value = prm.get_double("constant viscosity value");
    first_order_viscosity_coef =
      prm.get_double("first order viscosity coefficient");
    use_low_order_viscosity_for_first_time_step =
      prm.get_bool("use low order viscosity for first time step");
    entropy_residual_coef = prm.get_double("entropy residual coefficient");
    entropy_jump_coef = prm.get_double("entropy jump coefficient");
    entropy_viscosity_smoothing = prm.get("entropy viscosity smoothing");
    entropy_viscosity_smoothing_weight =
      prm.get_integer("entropy viscosity smoothing weight");
  }
  prm.leave_subsection();

  // output parameters
  prm.enter_subsection("output");
  {
    output_period = prm.get_integer("output period");
    max_transient_output_size = prm.get_double("max transient output size");
    output_mesh = prm.get_bool("output mesh");
    output_mass_matrix = prm.get_bool("output mass matrix");
    output_viscosity = prm.get_bool("output viscosity");
    output_viscosity_transient = prm.get_bool("output viscosity transient");
    output_exact_solution = prm.get_bool("output exact solution");
    exact_solution_refinement_level =
      prm.get_integer("exact solution refinement level");
    save_convergence_results = prm.get_bool("save convergence results");
    print_final_solution = prm.get_bool("print final solution");
    verbosity_level = prm.get_integer("verbosity level");
  }
  prm.leave_subsection();

  // FCT
  prm.enter_subsection("fct");
  {
    // antidiffusion
    std::string antidiffusion_string = prm.get("antidiffusion");
    if (antidiffusion_string == "limited")
      antidiffusion_type = AntidiffusionType::limited;
    else if (antidiffusion_string == "full")
      antidiffusion_type = AntidiffusionType::full;
    else if (antidiffusion_string == "none")
      antidiffusion_type = AntidiffusionType::none;
    else
      Assert(false, ExcNotImplemented());

    // synchronization
    std::string synchronization_string = prm.get("fct synchronization");
    if (synchronization_string == "none")
      fct_synchronization_type = FCTSynchronizationType::none;
    else if (synchronization_string == "min")
      fct_synchronization_type = FCTSynchronizationType::min;
    else if (synchronization_string == "compound")
      fct_synchronization_type = FCTSynchronizationType::compound;
    else
      Assert(false, ExcNotImplemented());

    // star states
    use_star_states_in_fct_bounds = prm.get_bool("use star states in fct bounds");
  }
  prm.leave_subsection();
}
