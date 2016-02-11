/** \file TransportParameters.h
 *  \brief Provides the header for the TransportParameters class.
 */
#ifndef TransportParameters_h
#define TransportParameters_h

#include <iostream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

#include "include/parameters/RunParameters.h"

using namespace dealii;

/** \class TransportParameters
 *  \brief Class for parameters related to a transport equation.
 */
template <int dim>
class TransportParameters : public RunParameters<dim>
{
public:
  TransportParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  std::string problem_name;
};

#include "src/parameters/TransportParameters.cc"

#endif
