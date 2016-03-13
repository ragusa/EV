/** \file EulerRunParameters.cc
 *  \brief Provides the function definitions for the EulerRunParameters class.
 */
using namespace dealii;

/** \fn EulerRunParameters::EulerRunParameters()
 *  \brief Constructor for the EulerRunParameters class.
 */
template <int dim>
EulerRunParameters::EulerRunParameters()
{
}

/**
 * \fn    EulerRunParameters::declare_euler_parameters(ParameterHandler
 * &parameter_handler)
 * \brief defines input parameters
 * \param parameter_handler parameter handler for the Euler class
 */
template <int dim>
void EulerRunParameters::declare_euler_parameters(
  ParameterHandler & parameter_handler)
{
  // declare conservation law parameters
  RunParameters::declare_run_parameters(parameter_handler);

  // problem
  parameter_handler.enter_subsection("problem");
  {
    parameter_handler.declare_entry("problem id",
                                    "0",
                                    Patterns::Integer(),
                                    "ID for description of the problem");
  }
  parameter_handler.leave_subsection();

  // numerical parameters
  parameter_handler.enter_subsection("numerics");
  {
    parameter_handler.declare_entry(
      "prandtl",
      "0.0",
      Patterns::Double(),
      "Prandtl number; used in relating viscosities nu and kappa");
  }
  parameter_handler.leave_subsection();
}

/**
 * \fn    EulerRunParameters::get_euler_parameters(ParameterHandler
 * &parameter_handler)
 * \brief get input parameters from parameter handler
 * \param parameter_handler parameter handler for the Euler class
 */
template <int dim>
void EulerRunParameters::get_euler_parameters(
  ParameterHandler & parameter_handler)

{
  // get conservation law parameters
  this->get_run_parameters(parameter_handler);
  this->n_components = n_euler_components;

  // problem
  parameter_handler.enter_subsection("problem");
  {
    problem_id = parameter_handler.get_integer("problem id");
  }
  parameter_handler.leave_subsection();

  // numerical parameters
  parameter_handler.enter_subsection("numerics");
  {
    prandtl = parameter_handler.get_double("prandtl");
  }
  parameter_handler.leave_subsection();
}
