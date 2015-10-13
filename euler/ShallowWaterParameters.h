/**
 * \file ShallowWaterParameters.h
 * \brief Provides the header for the ShallowWaterParameters class.
 */
#ifndef ShallowWaterParameters_h
#define ShallowWaterParameters_h

#include <iostream>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

using namespace dealii;

/**
 * \class ShallowWaterParameters
 * \brief Class for parameters related to the shallow water equations.
 */
template <int dim>
class ShallowWaterParameters : public ConservationLawParameters<dim>
{
public:
  ShallowWaterParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  unsigned int problem_id;

  bool use_local_entropy_normalization;

  bool multiply_low_order_viscosity_by_froude;
};

#include "ShallowWaterParameters.cc"

#endif
