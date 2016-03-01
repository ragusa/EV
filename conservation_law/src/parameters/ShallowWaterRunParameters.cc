/**
 * \file ShallowWaterRunParameters.cc
 * \brief Provides the function definitions for the ShallowWaterRunParameters class.
 */
using namespace dealii;

/**
 * \brief Constructor for the ShallowWaterRunParameters class
 */
template <int dim>
ShallowWaterRunParameters<dim>::ShallowWaterRunParameters()
{
}

/**
 * \brief Defines input parameters.
 *
 * \param[out] parameter_handler parameter handler for the ShallowWater class
 */
template <int dim>
void ShallowWaterRunParameters<dim>::declare_parameters(
  ParameterHandler & parameter_handler)
{
  // declare conservation law parameters
  RunParameters<dim>::declare_run_parameters(parameter_handler);

  // problem
  parameter_handler.enter_subsection("problem");
  {
    parameter_handler.declare_entry(
      "problem name", "default", Patterns::Anything(), "Problem name");
  }
  parameter_handler.leave_subsection();

  // artificial viscosity
  parameter_handler.enter_subsection("artificial viscosity");
  {
    parameter_handler.declare_entry("entropy normalization",
                                    "average",
                                    Patterns::Anything(),
                                    "Option for normalization of entropy"
                                    " viscosity");
    parameter_handler.declare_entry("constant entropy normalization coefficient",
                                    "1.0",
                                    Patterns::Double(),
                                    "Entropy viscosity normalization value if"
                                    " constant normalization option is chosen");
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
void ShallowWaterRunParameters<dim>::get_parameters(
  ParameterHandler & parameter_handler)
{
  // get conservation law parameters
  this->get_run_parameters(parameter_handler);

  // problem
  parameter_handler.enter_subsection("problem");
  {
    problem_name = parameter_handler.get("problem name");
  }
  parameter_handler.leave_subsection();

  // artificial viscosity
  parameter_handler.enter_subsection("artificial viscosity");
  {
    entropy_normalization = parameter_handler.get("entropy normalization");
    constant_entropy_normalization_coefficient =
      parameter_handler.get_double("constant entropy normalization coefficient");
    multiply_low_order_viscosity_by_froude =
      parameter_handler.get_bool("multiply low order viscosity by froude");
  }
  parameter_handler.leave_subsection();
}
