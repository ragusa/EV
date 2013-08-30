/** \file EulerEquationsParameters.h
 *  \brief Provides the header for the EulerEquationsParameters class.
 */
#ifndef EulerEquationsParameters_h
#define EulerEquationsParameters_h

#include <iostream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

using namespace dealii;

/** \class EulerEquationsParameters
 *  \brief Defines and retrieves input parameters related to Euler equations
 */
template <int dim>
class EulerEquationsParameters
{
  public:
    EulerEquationsParameters();

    static void declare_parameters (ParameterHandler &parameter_handler);
           void get_parameters     (ParameterHandler &parameter_handler);

    static const int n_components = dim + 2;

    double input1;
    double input2;
    double input3;

    std::vector<std::string> initial_conditions_expressions;
};

#include "EulerEquationsParameters.cc"

#endif
