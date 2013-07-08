#include <deal.II/base/parameter_handler.h>

#include "EulerEquationsBaseParameters.h"

using namespace dealii;

/**
 * \fn    EulerEquationsBaseParameters::declare_parameters
 * \brief defines input parameters
 */
template<int dim>
void EulerEquationsBaseParameters<dim>::declare_parameters(
  ParameterHandler &parameter_handler)
{
  parameter_handler.declare_entry("First input parameter", "5.0",
    Patterns::Double(),"Description of first input parameter");
  parameter_handler.declare_entry("Second input parameter", "5.0",
    Patterns::Double(),"Description of second input parameter");
  parameter_handler.declare_entry("Third input parameter", "5.0",
    Patterns::Double(),"Description of third input parameter");
}

/**
 * \fn    EulerEquationsBaseParameters::get_parameters
 * \brief get input parameters from parameter handler
 */
template<int dim>
void EulerEquationsBaseParameters<dim>::get_parameters(
  ParameterHandler &parameter_handler)
{
  input1 = parameter_handler.get_double("First input parameter");
  input2 = parameter_handler.get_double("Second input parameter");
  input3 = parameter_handler.get_double("Third input parameter");
}
