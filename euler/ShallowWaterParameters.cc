/**
 * \file ShallowWaterParameters.cc
 * \brief Provides the function definitions for the ShallowWaterParameters class.
 */
using namespace dealii;

/**
 * \brief Constructor for the ShallowWaterParameters class
 */
template <int dim>
ShallowWaterParameters<dim>::ShallowWaterParameters()
{
}

/**
 * \brief Defines input parameters.
 *
 * \param[out] parameter_handler parameter handler for the ShallowWater class
 */
template <int dim>
void ShallowWaterParameters<dim>::declare_parameters(
  ParameterHandler & parameter_handler)
{
  // declare conservation law parameters
  ConservationLawParameters<dim>::declare_conservation_law_parameters(
    parameter_handler);

  // problem
  parameter_handler.enter_subsection("problem");
  {
    parameter_handler.declare_entry("problem id",
                                    "0",
                                    Patterns::Integer(),
                                    "ID for description of the problem");
  }
  parameter_handler.leave_subsection();

  // artificial viscosity
  parameter_handler.enter_subsection("artificial viscosity");
  {
    parameter_handler.declare_entry("use local entropy normalization",
                                    "false",
                                    Patterns::Bool(),
                                    "Option to use a local normalization for"
                                    " the entropy residual");
    parameter_handler.declare_entry("multiply low order viscosity by froude",
                                    "false",
                                    Patterns::Bool(),
                                    "Option to multiply the low-order viscisity"
                                    " by the local Froude number");
  }
  parameter_handler.leave_subsection();
}

/**
 * \brief Gets input parameters from parameter handler.
 *
 * \param[in] parameter_handler parameter handler for the ShallowWater class
 */
template <int dim>
void ShallowWaterParameters<dim>::get_parameters(
  ParameterHandler & parameter_handler)
{
  // get conservation law parameters
  this->get_conservation_law_parameters(parameter_handler);
  this->n_components = dim + 1;

  // problem
  parameter_handler.enter_subsection("problem");
  {
    problem_id = parameter_handler.get_integer("problem id");
  }
  parameter_handler.leave_subsection();

  // artificial viscosity
  parameter_handler.enter_subsection("artificial viscosity");
  {
    use_local_entropy_normalization =
      parameter_handler.get_bool("use local entropy normalization");
    multiply_low_order_viscosity_by_froude =
      parameter_handler.get_bool("multiply low order viscosity by froude");
  }
  parameter_handler.leave_subsection();
}
