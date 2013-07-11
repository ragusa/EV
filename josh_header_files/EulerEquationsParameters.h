/*
 EulerEquationsParameters class

 contains input data relating to Euler equations
*/

#ifndef EulerEquationsParameters_h
#define EulerEquationsParameters_h

#include <iostream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

using namespace dealii;

template <int dim>
class EulerEquationsParameters
{
  public:
    EulerEquationsParameters();

    void declare_parameters(ParameterHandler &parameter_handler);
    void get_parameters    (ParameterHandler &parameter_handler);

    static const int n_components = dim + 2;

    double input1;
    double input2;
    double input3;
};

/* The source file must be included here because the class is
   a template; the compiler must be able to see the class template
   declaration AND definition.
*/
#include "EulerEquationsParameters.cc"

#endif
