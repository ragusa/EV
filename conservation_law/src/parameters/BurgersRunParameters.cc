/** \file BurgersRunParameters.cc
 *  \brief Provides the function definitions for the BurgersRunParameters class.
 */
using namespace dealii;

/**
 * \brief Constructor for the BurgersRunParameters class
 */
template <int dim>
BurgersRunParameters::BurgersRunParameters()
{
}

/**
 * \brief defines input parameters
 * \param parameter_handler parameter handler for the Burgers class
 */
template <int dim>
void BurgersRunParameters::declare_parameters(
  ParameterHandler & parameter_handler)
{
  // declare conservation law parameters
  RunParameters::declare_run_parameters(parameter_handler);

  // problem
  parameter_handler.enter_subsection("problem");
  {
    parameter_handler.declare_entry(
      "problem name", "default", Patterns::Anything(), "Problem name");
  }
  parameter_handler.leave_subsection();
}

/**
 * \brief gets input parameters from parameter handler
 * \param parameter_handler parameter handler for the Burgers class
 */
template <int dim>
void BurgersRunParameters::get_parameters(ParameterHandler & parameter_handler)
{
  // get conservation law parameters
  this->get_run_parameters(parameter_handler);
  this->n_components = 1;

  // problem
  parameter_handler.enter_subsection("problem");
  {
    problem_name = parameter_handler.get("problem name");
  }
  parameter_handler.leave_subsection();
}
