/** \file EulerParameters.h
 *  \brief Provides the header for the EulerParameters class.
 */
#ifndef EulerParameters_h
#define EulerParameters_h

#include <iostream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

using namespace dealii;

/** \class EulerParameters
 *  \brief Declares and retrieves parameters related to
 *         the Euler equations.
 */
template <int dim>
class EulerParameters : public ConservationLawParameters<dim>
{
  public:
    EulerParameters();

    static void declare_euler_parameters (ParameterHandler &parameter_handler);
           void get_euler_parameters     (ParameterHandler &parameter_handler);

    static const int n_euler_components = 1;

    unsigned int problem_id;
};

#include "EulerParameters.cc"

#endif
