/**
 * \file TransportProblemParameters.h
 * \brief Provides the header for the TransportProblemParameters class.
 */
#ifndef TransportProblemParameters_h
#define TransportProblemParameters_h

#include <deal.II/base/parameter_handler.h>
#include "include/parameters/ProblemParameters.h"

using namespace dealii;

/**
 * \class TransportProblemParameters
 * \brief Class for parameters related to transport problems.
 */
template <int dim>
class TransportProblemParameters : public ProblemParameters<dim>
{
public:
  TransportProblemParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  std::string dirichlet_function;
  std::string initial_condition;
  std::string exact_solution;
};

#include "src/parameters/TransportProblemParameters.cc"

#endif
