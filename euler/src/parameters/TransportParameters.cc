/** \file TransportParameters.cc
 *  \brief Provides the function definitions for the TransportParameters class.
 */
using namespace dealii;

/**
 * \brief Constructor for the TransportParameters class
 */
template <int dim>
TransportParameters<dim>::TransportParameters()
{
}

/**
 * \brief defines input parameters
 * \param parameter_handler parameter handler for the Transport class
 */
template <int dim>
void TransportParameters<dim>::declare_parameters(
  ParameterHandler & parameter_handler)
{
  // declare conservation law parameters
  ConservationLawParameters<dim>::declare_conservation_law_parameters(
    parameter_handler);

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
 * \param parameter_handler parameter handler for the Transport class
 */
template <int dim>
void TransportParameters<dim>::get_parameters(ParameterHandler & parameter_handler)
{
  // get conservation law parameters
  this->get_conservation_law_parameters(parameter_handler);
  this->n_components = 1;

  // problem
  parameter_handler.enter_subsection("problem");
  {
    problem_name = parameter_handler.get("problem name");
  }
  parameter_handler.leave_subsection();
}
