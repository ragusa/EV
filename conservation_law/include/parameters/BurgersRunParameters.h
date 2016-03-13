/** \file BurgersRunParameters.h
 *  \brief Provides the header for the BurgersRunParameters class.
 */
#ifndef BurgersRunParameters_h
#define BurgersRunParameters_h

#include <iostream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

#include "include/parameters/RunParameters.h"

using namespace dealii;

/** \class BurgersRunParameters
 *  \brief Class for parameters related to the Burgers equation.
 */
template <int dim>
class BurgersRunParameters : public RunParameters
{
public:
  BurgersRunParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  std::string problem_name;
};

#include "src/parameters/BurgersRunParameters.cc"

#endif
