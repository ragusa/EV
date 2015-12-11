/**
 * \file TransportProblemParameters.cc
 * \brief Provides the function definitions for the TransportProblemParameters
 * class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
TransportProblemParameters<dim>::TransportProblemParameters()
{
}

/**
 * \brief Declares parameters
 *
 * \param[out] parameter_handler parameter handler for the Transport class
 */
template <int dim>
void TransportProblemParameters<dim>::declare_parameters(
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
      "initial condition", "0", Patterns::Anything(), "Initial conditions");
  }
  parameter_handler.leave_subsection();

  // exact solution
  parameter_handler.enter_subsection("exact solution");
  {
    parameter_handler.declare_entry(
      "exact solution", "1", Patterns::Anything(), "Exact solution");
  }
  parameter_handler.leave_subsection();
}

/**
 * \brief Gets parameters from parameter handler.
 *
 * \param[in] parameter_handler parameter handler
 */
template <int dim>
void TransportProblemParameters<dim>::get_parameters(
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
    initial_conditions = parameter_handler.get("initial condition");
  }
  parameter_handler.leave_subsection();

  // exact solution
  parameter_handler.enter_subsection("exact solution");
  {
    exact_solution = parameter_handler.get("exact solution");
  }
  parameter_handler.leave_subsection();
}
