/**
 * \file TransportRunParameters.cc
 * \brief Provides the function definitions for the TransportRunParameters class.
 */
using namespace dealii;

/**
 * \brief Constructor.
 */
// template <int dim>
// TransportRunParameters::TransportRunParameters()
TransportRunParameters::TransportRunParameters() {}
/**
 * \brief Defines input parameters.
 *
 * \param[in] parameter_handler parameter handler for the Transport class
 */
// template <int dim>
// void TransportRunParameters::declare_parameters(
void TransportRunParameters::declare_parameters(
  ParameterHandler & parameter_handler)
{
  // declare base run parameters
  // RunParameters::declare_run_parameters(parameter_handler);
  RunParameters::declare_run_parameters(parameter_handler);
}

/**
 * \brief Gets input parameters from parameter handler.
 *
 * \param[in] parameter_handler parameter handler for the Transport class
 */
// template <int dim>
// void TransportRunParameters::get_parameters(
void TransportRunParameters::get_parameters(ParameterHandler & parameter_handler)
{
  // get conservation law parameters
  this->get_run_parameters(parameter_handler);
}
