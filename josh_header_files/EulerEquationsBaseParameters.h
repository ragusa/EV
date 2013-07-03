/*
 EulerEquationsBaseParameters class

 contains input data relating to Euler equations
*/

#ifndef EulerEquationsBaseParameters_h
#define EulerEquationsBaseParameters_h

#include <deal.II/base/parameter_handler.h>

class EulerEquationsBaseParameters
{
  public:
    EulerEquationsBaseParameters();
    static void declare_parameters(ParameterHandler &parameter_handler);
    void        get_parameters    (ParameterHandler &parameter_handler);

    double input1;
    double input2;
    double input3;
};

#endif
