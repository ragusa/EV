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
    prm.declare_entry("diffusion type",
                      "laplacian",
                      Patterns::Anything(),
                      "choice for type of diffusion");
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
      time_step_size_method = TimeStepSizeMethod::constant_dt;
    else if (time_choice == "cfl_condition")
      time_step_size_method = TimeStepSizeMethod::cfl_condition;
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
    if (rk_choice == "ERK1")
      time_discretization = TemporalDiscretization::ERK1;
    else if (rk_choice == "ERK2")
      time_discretization = TemporalDiscretization::ERK2;
    else if (rk_choice == "ERK3")
      time_discretization = TemporalDiscretization::ERK3;
    else if (rk_choice == "ERK4")
      time_discretization = TemporalDiscretization::ERK4;
    else if (rk_choice == "SDIRK22")
      time_discretization = TemporalDiscretization::SDIRK22;
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
    const std::string viscosity_choice = prm.get("viscosity type");
    if (viscosity_choice == "none")
      viscosity_type = ViscosityType::none;
    else if (viscosity_choice == "constant")
      viscosity_type = ViscosityType::constant;
    else if (viscosity_choice == "low")
      viscosity_type = ViscosityType::low;
    else if (viscosity_choice == "DMP_low")
      viscosity_type = ViscosityType::DMP_low;
    else if (viscosity_choice == "DI_low")
      viscosity_type = ViscosityType::DI_low;
    else if (viscosity_choice == "entropy")
      viscosity_type = ViscosityType::entropy;
    else
      Assert(false, ExcNotImplemented());

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

    const std::string diffusion_choice = prm.get("diffusion type");
    if (diffusion_choice == "none")
      diffusion_type = DiffusionType::none;
    else if (diffusion_choice == "laplacian")
      diffusion_type = DiffusionType::laplacian;
    else if (diffusion_choice == "graphtheoretic")
      diffusion_type = DiffusionType::graphtheoretic;
    else
      ExcNotImplemented();
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
}