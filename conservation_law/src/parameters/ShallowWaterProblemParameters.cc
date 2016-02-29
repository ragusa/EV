/**
 * \file ShallowWaterProblemParameters.cc
 * \brief Provides the function definitions for the ShallowWaterProblemParameters
 *        class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ShallowWaterProblemParameters<dim>::ShallowWaterProblemParameters()
{
}

/**
 * \brief Declares parameters
 *
 * \param[out] parameter_handler parameter handler for the ShallowWater class
 */
template <int dim>
void ShallowWaterProblemParameters<dim>::declare_parameters(
  ParameterHandler & parameter_handler)
{
  // declare common problem parameters
  ProblemParameters<dim>::declare_common_parameters(parameter_handler);

  // bathymetry
  parameter_handler.enter_subsection("bathymetry");
  {
    parameter_handler.declare_entry(
      "bathymetry function", "0", Patterns::Anything(), "Bathymetry function");
  }
  parameter_handler.leave_subsection();

  // boundary conditions
  parameter_handler.enter_subsection("boundary conditions");
  {
    parameter_handler.declare_entry("dirichlet function height",
                                    "0",
                                    Patterns::Anything(),
                                    "Dirichlet function for height");
    parameter_handler.declare_entry("dirichlet function momentumx",
                                    "0",
                                    Patterns::Anything(),
                                    "Dirichlet function for x-momentum");
    parameter_handler.declare_entry("dirichlet function momentumy",
                                    "0",
                                    Patterns::Anything(),
                                    "Dirichlet function for y-momentum");
  }
  parameter_handler.leave_subsection();

  // initial conditions
  parameter_handler.enter_subsection("initial conditions");
  {
    parameter_handler.declare_entry("initial conditions height",
                                    "0",
                                    Patterns::Anything(),
                                    "Initial conditions for height");
    parameter_handler.declare_entry("initial conditions momentumx",
                                    "0",
                                    Patterns::Anything(),
                                    "Initial conditions for x-momentum");
    parameter_handler.declare_entry("initial conditions momentumy",
                                    "0",
                                    Patterns::Anything(),
                                    "Initial conditions for y-momentum");
  }
  parameter_handler.leave_subsection();

  // exact solution
  parameter_handler.enter_subsection("exact solution");
  {
    parameter_handler.declare_entry("exact solution height",
                                    "1",
                                    Patterns::Anything(),
                                    "Exact solution for height");
    parameter_handler.declare_entry("exact solution momentumx",
                                    "1",
                                    Patterns::Anything(),
                                    "Exact solution for x-momentum");
    parameter_handler.declare_entry("exact solution momentumy",
                                    "1",
                                    Patterns::Anything(),
                                    "Exact solution for y-momentum");
  }
  parameter_handler.leave_subsection();

  // constants
  parameter_handler.enter_subsection("constants");
  {
    parameter_handler.declare_entry(
      "gravity", "1.0", Patterns::Double(), "Acceleration due to gravity");
    parameter_handler.declare_entry(
      "x_interface", "0.0", Patterns::Double(), "x-position of interface");
    parameter_handler.declare_entry(
      "h_left", "1.0", Patterns::Double(), "Left height value");
    parameter_handler.declare_entry(
      "h_right", "1.0", Patterns::Double(), "Right height value");
    parameter_handler.declare_entry(
      "h_unperturbed", "1.0", Patterns::Double(), "Unperturbed height value");
    parameter_handler.declare_entry(
      "h_perturbed", "1.0", Patterns::Double(), "Perturbed height value");
    parameter_handler.declare_entry(
      "u_left", "0.0", Patterns::Double(), "Left x-velocity value");
    parameter_handler.declare_entry(
      "u_right", "0.0", Patterns::Double(), "Right x-velocity value");
    parameter_handler.declare_entry("bump_x_center",
                                    "0.0",
                                    Patterns::Double(),
                                    "Center of bump in x-direction");
    parameter_handler.declare_entry("bump_y_center",
                                    "0.0",
                                    Patterns::Double(),
                                    "Center of bump in y-direction");
    parameter_handler.declare_entry(
      "bump_x_width", "0.0", Patterns::Double(), "Width of bump in x-direction");
    parameter_handler.declare_entry(
      "bump_y_width", "0.0", Patterns::Double(), "Width of bump in y-direction");
    parameter_handler.declare_entry(
      "bump_height", "0.0", Patterns::Double(), "Height of bump");
    parameter_handler.declare_entry("perturbation_x_center",
                                    "0.0",
                                    Patterns::Double(),
                                    "Center of perturbation in x-direction");
    parameter_handler.declare_entry("perturbation_y_center",
                                    "0.0",
                                    Patterns::Double(),
                                    "Center of perturbation in y-direction");
    parameter_handler.declare_entry("perturbation_x_width",
                                    "0.0",
                                    Patterns::Double(),
                                    "Width of perturbation in x-direction");
    parameter_handler.declare_entry("perturbation_y_width",
                                    "0.0",
                                    Patterns::Double(),
                                    "Width of perturbation in y-direction");
  }
  parameter_handler.leave_subsection();
}

/**
 * \brief Gets parameters from parameter handler.
 *
 * \param[in] parameter_handler parameter handler
 */
template <int dim>
void ShallowWaterProblemParameters<dim>::get_parameters(
  ParameterHandler & parameter_handler)
{
  // get common problem parameters
  this->get_common_parameters(parameter_handler);

  // bathymetry
  parameter_handler.enter_subsection("bathymetry");
  {
    bathymetry_function = parameter_handler.get("bathymetry function");
  }
  parameter_handler.leave_subsection();

  // boundary conditions
  parameter_handler.enter_subsection("boundary conditions");
  {
    dirichlet_function_height =
      parameter_handler.get("dirichlet function height");
    dirichlet_function_momentumx =
      parameter_handler.get("dirichlet function momentumx");
    dirichlet_function_momentumy =
      parameter_handler.get("dirichlet function momentumy");
  }
  parameter_handler.leave_subsection();

  // initial conditions
  parameter_handler.enter_subsection("initial conditions");
  {
    initial_conditions_height =
      parameter_handler.get("initial conditions height");
    initial_conditions_momentumx =
      parameter_handler.get("initial conditions momentumx");
    initial_conditions_momentumy =
      parameter_handler.get("initial conditions momentumy");
  }
  parameter_handler.leave_subsection();

  // exact solution
  parameter_handler.enter_subsection("exact solution");
  {
    exact_solution_height = parameter_handler.get("exact solution height");
    exact_solution_momentumx = parameter_handler.get("exact solution momentumx");
    exact_solution_momentumy = parameter_handler.get("exact solution momentumy");
  }
  parameter_handler.leave_subsection();

  // constants
  parameter_handler.enter_subsection("constants");
  {
    gravity = parameter_handler.get_double("gravity");
    x_interface = parameter_handler.get_double("x_interface");
    h_left = parameter_handler.get_double("h_left");
    h_right = parameter_handler.get_double("h_right");
    h_unperturbed = parameter_handler.get_double("h_unperturbed");
    h_perturbed = parameter_handler.get_double("h_perturbed");
    u_left = parameter_handler.get_double("u_left");
    u_right = parameter_handler.get_double("u_right");
    bump_x_center = parameter_handler.get_double("bump_x_center");
    bump_y_center = parameter_handler.get_double("bump_y_center");
    bump_x_width = parameter_handler.get_double("bump_x_width");
    bump_y_width = parameter_handler.get_double("bump_y_width");
    bump_height = parameter_handler.get_double("bump_height");
    perturbation_x_center = parameter_handler.get_double("perturbation_x_center");
    perturbation_y_center = parameter_handler.get_double("perturbation_y_center");
    perturbation_x_width = parameter_handler.get_double("perturbation_x_width");
    perturbation_y_width = parameter_handler.get_double("perturbation_y_width");
  }
  parameter_handler.leave_subsection();
}
