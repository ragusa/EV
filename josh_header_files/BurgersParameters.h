/** \file BurgersParameters.h
 *  \brief Provides the header for the BurgersParameters class.
 */
#ifndef BurgersParameters_h
#define BurgersParameters_h

#include <iostream>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

using namespace dealii;

/** \class BurgersParameters
 *  \brief Declares and retrieves parameters related to
 *         Burgers equation.
 */
template <int dim>
class BurgersParameters : public ConservationLawParameters<dim>
{
  public:
    BurgersParameters();

    static void declare_burgers_parameters (ParameterHandler &parameter_handler);
           void get_burgers_parameters     (ParameterHandler &parameter_handler);

    static const int n_burgers_components = 1;

    unsigned int problem_id;
};

#include "BurgersParameters.cc"

#endif
