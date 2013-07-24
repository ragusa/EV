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
class BurgersParameters
{
  public:
    BurgersParameters();

    static void declare_parameters (ParameterHandler &parameter_handler);
           void get_parameters     (ParameterHandler &parameter_handler);

    static const int n_components = 1;

    std::vector<std::string> initial_conditions_expressions;
};

#include "BurgersParameters.cc"

#endif
