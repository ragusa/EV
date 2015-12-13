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

  // physics
  parameter_handler.enter_subsection("physics");
  {
    parameter_handler.declare_entry(
      "transport speed", "1.0", Patterns::Double(), "Transport speed");
    parameter_handler.declare_entry("transport direction x",
                                    "1.0",
                                    Patterns::Double(),
                                    "x-component of transport direction");
    parameter_handler.declare_entry("transport direction y",
                                    "0.0",
                                    Patterns::Double(),
                                    "y-component of transport direction");
    parameter_handler.declare_entry("transport direction z",
                                    "0.0",
                                    Patterns::Double(),
                                    "z-component of transport direction");
    parameter_handler.declare_entry(
      "cross section", "0", Patterns::Anything(), "Cross section");
    parameter_handler.declare_entry(
      "source", "0", Patterns::Anything(), "Source");
  }
  parameter_handler.leave_subsection();

  // constants
  parameter_handler.enter_subsection("constants");
  {
    parameter_handler.declare_entry("incoming value",
                                    "1.0",
                                    Patterns::Double(),
                                    "incoming value for function parsers");
    parameter_handler.declare_entry("cross section value",
                                    "0.0",
                                    Patterns::Double(),
                                    "cross section for function parsers");
    parameter_handler.declare_entry(
      "source value", "0.0", Patterns::Double(), "source for function parsers");
  }
  parameter_handler.leave_subsection();

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

  // physics
  parameter_handler.enter_subsection("physics");
  {
    transport_speed = parameter_handler.get_double("transport speed");
    transport_direction_x = parameter_handler.get_double("transport direction x");
    transport_direction_y = parameter_handler.get_double("transport direction y");
    transport_direction_z = parameter_handler.get_double("transport direction z");
    cross_section_string = parameter_handler.get("cross section");
    source_string = parameter_handler.get("source");
  }
  parameter_handler.leave_subsection();

  // constants
  parameter_handler.enter_subsection("constants");
  {
    incoming_value = parameter_handler.get_double("incoming value");
    cross_section_value = parameter_handler.get_double("cross section value");
    source_value = parameter_handler.get_double("source value");
  }
  parameter_handler.leave_subsection();

  // boundary conditions
  parameter_handler.enter_subsection("boundary conditions");
  {
    dirichlet_function = parameter_handler.get("dirichlet function");
  }
  parameter_handler.leave_subsection();

  // initial conditions
  parameter_handler.enter_subsection("initial conditions");
  {
    initial_condition = parameter_handler.get("initial condition");
  }
  parameter_handler.leave_subsection();

  // exact solution
  parameter_handler.enter_subsection("exact solution");
  {
    exact_solution = parameter_handler.get("exact solution");
  }
  parameter_handler.leave_subsection();
}
