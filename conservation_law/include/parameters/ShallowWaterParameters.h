/**
 * \file ShallowWaterParameters.h
 * \brief Provides the header for the ShallowWaterParameters class.
 */
#ifndef ShallowWaterParameters_h
#define ShallowWaterParameters_h

#include <iostream>
#include <string>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

#include "include/parameters/RunParameters.h"

using namespace dealii;

/**
 * \class ShallowWaterParameters
 * \brief Class for parameters related to the shallow water equations.
 */
template <int dim>
class ShallowWaterParameters : public RunParameters<dim>
{
public:
  ShallowWaterParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  std::string problem_name;

  std::string entropy_normalization;

  double constant_entropy_normalization_coefficient;

  bool multiply_low_order_viscosity_by_froude;
};

#include "src/parameters/ShallowWaterParameters.cc"

#endif
