/*
 EulerEquationsBaseParameters class

 contains input data relating to Euler equations
*/

#ifndef EulerEquationsBaseParameters_h
#define EulerEquationsBaseParameters_h

#include <iostream>

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

template <int dim>
class EulerEquationsBaseParameters
{
  public:
    void declare_parameters(ParameterHandler &parameter_handler);
    void get_parameters    (ParameterHandler &parameter_handler);

    double input1;
    double input2;
    double input3;
};

/* The source file must be included here because the class is
   a template; the compiler must be able to see the class template
   declaration AND definition.
*/
#include "EulerEquationsBaseParameters.cc"

#endif
