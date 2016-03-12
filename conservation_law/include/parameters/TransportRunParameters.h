/** \file TransportRunParameters.h
 *  \brief Provides the header for the TransportRunParameters class.
 */
#ifndef TransportRunParameters_h
#define TransportRunParameters_h

#include <iostream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

#include "include/parameters/RunParameters.h"

using namespace dealii;

/** \class TransportRunParameters
 *  \brief Class for parameters related to a transport equation.
 */
//template <int dim>
//class TransportRunParameters : public RunParameters<dim>
class TransportRunParameters : public RunParameters
{
public:
  TransportRunParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);
};

#include "src/parameters/TransportRunParameters.cc"

#endif
