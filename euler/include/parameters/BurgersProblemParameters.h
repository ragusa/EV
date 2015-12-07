/**
 * \file BurgersProblemParameters.h
 * \brief Provides the header for the BurgersProblemParameters class.
 */
#ifndef BurgersProblemParameters_h
#define BurgersProblemParameters_h

using namespace dealii;

/**
 * \class BurgersProblemParameters
 * \brief Class for parameters related to Burgers problems.
 */
template <int dim>
class BurgersProblemParameters
{
public:
  BurgersProblemParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  std::string dirichlet_function;
  std::string initial_conditions;
  std::string exact_solution;

  double x_interface;
  double u_left;
  double u_right;
};

#include "src/parameters/BurgersProblemParameters.cc"

#endif
