/** \file EulerRunParameters.h
 *  \brief Provides the header for the EulerRunParameters class.
 */
#ifndef EulerRunParameters_h
#define EulerRunParameters_h

#include <iostream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

#include "include/parameters/RunParameters.h"

using namespace dealii;

/** \class EulerRunParameters
 *  \brief Class for parameters related to the Euler equations.
 */
template <int dim>
class EulerRunParameters : public RunParameters
{
public:
  EulerRunParameters();

  static void declare_euler_parameters(ParameterHandler & parameter_handler);

  void get_euler_parameters(ParameterHandler & parameter_handler);

  static const int n_euler_components = dim + 2;

  unsigned int problem_id;

  double prandtl; // Prandtl number, used in relation of viscosities nu and kappa
};

#include "src/parameters/EulerRunParameters.cc"

#endif
