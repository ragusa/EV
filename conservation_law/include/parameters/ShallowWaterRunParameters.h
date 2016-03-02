/**
 * \file ShallowWaterRunParameters.h
 * \brief Provides the header for the ShallowWaterRunParameters class.
 */
#ifndef ShallowWaterRunParameters_h
#define ShallowWaterRunParameters_h

#include <iostream>
#include <string>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

#include "include/parameters/RunParameters.h"

using namespace dealii;

/**
 * \class ShallowWaterRunParameters
 * \brief Class for parameters related to the shallow water equations.
 */
template <int dim>
class ShallowWaterRunParameters : public RunParameters<dim>
{
public:
  ShallowWaterRunParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  std::string entropy_normalization;

  double constant_entropy_normalization_coefficient;

  bool multiply_low_order_viscosity_by_froude;
};

#include "src/parameters/ShallowWaterRunParameters.cc"

#endif
