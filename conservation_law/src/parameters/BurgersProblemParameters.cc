/**
 * \file BurgersProblemParameters.cc
 * \brief Provides the function definitions for the BurgersProblemParameters
 * class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
BurgersProblemParameters<dim>::BurgersProblemParameters()
{
}

/**
 * \brief Declares parameters
 *
 * \param[out] parameter_handler parameter handler for the Burgers class
 */
template <int dim>
void BurgersProblemParameters<dim>::declare_parameters(
  ParameterHandler & parameter_handler)
{
  // declare common problem parameters
  ProblemParameters<dim>::declare_common_parameters(parameter_handler);

  // boundary conditions
  parameter_handler.enter_subsection("boundary conditions");
  {
    parameter_handler.declare_entry(
      "dirichlet function", "0", Patterns::Anything(), "Dirichlet function");
  }
  parameter_handler.leave_subsection();

  // initial conditions
  parameter_handler.enter_subsection("initial conditions");
  {
    parameter_handler.declare_entry(
      "initial conditions", "0", Patterns::Anything(), "Initial conditions");
  }
  parameter_handler.leave_subsection();

  // exact solution
  parameter_handler.enter_subsection("exact solution");
  {
    parameter_handler.declare_entry(
      "exact solution", "1", Patterns::Anything(), "Exact solution");
  }
  parameter_handler.leave_subsection();

  // constants
  parameter_handler.enter_subsection("constants");
  {
    parameter_handler.declare_entry(
      "x_interface", "0.0", Patterns::Double(), "x-position of interface");
    parameter_handler.declare_entry(
      "u_left", "0.0", Patterns::Double(), "Left x-velocity value");
    parameter_handler.declare_entry(
      "u_right", "0.0", Patterns::Double(), "Right x-velocity value");
  }
  parameter_handler.leave_subsection();
}

/**
 * \brief Gets parameters from parameter handler.
 *
 * \param[in] parameter_handler parameter handler
 */
template <int dim>
void BurgersProblemParameters<dim>::get_parameters(
  ParameterHandler & parameter_handler)
{
  // get common problem parameters
  this->get_common_parameters(parameter_handler);

  // boundary conditions
  parameter_handler.enter_subsection("boundary conditions");
  {
    dirichlet_function = parameter_handler.get("dirichlet function");
  }
  parameter_handler.leave_subsection();

  // initial conditions
  parameter_handler.enter_subsection("initial conditions");
  {
    initial_conditions = parameter_handler.get("initial conditions");
  }
  parameter_handler.leave_subsection();

  // exact solution
  parameter_handler.enter_subsection("exact solution");
  {
    exact_solution = parameter_handler.get("exact solution");
  }
  parameter_handler.leave_subsection();

  // constants
  parameter_handler.enter_subsection("constants");
  {
    x_interface = parameter_handler.get_double("x_interface");
    u_left = parameter_handler.get_double("u_left");
    u_right = parameter_handler.get_double("u_right");
  }
  parameter_handler.leave_subsection();
}
