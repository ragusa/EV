/**
 * \file ShallowWaterProblemParameters.h
 * \brief Provides the header for the ShallowWaterProblemParameters class.
 */
#ifndef ShallowWaterProblemParameters_h
#define ShallowWaterProblemParameters_h

#include <iostream>
#include <string>
#include <deal.II/base/parameter_handler.h>
#include "include/parameters/ProblemParameters.h"

using namespace dealii;

/**
 * \class ShallowWaterProblemParameters
 * \brief Class for parameters related to the shallow water problems.
 */
template <int dim>
class ShallowWaterProblemParameters : public ProblemParameters<dim>
{
public:
  ShallowWaterProblemParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  std::string bathymetry_function;

  std::string dirichlet_function_height;
  std::string dirichlet_function_momentumx;
  std::string dirichlet_function_momentumy;

  std::string initial_conditions_height;
  std::string initial_conditions_momentumx;
  std::string initial_conditions_momentumy;

  std::string exact_solution_height;
  std::string exact_solution_momentumx;
  std::string exact_solution_momentumy;

  double gravity;
  double x_interface;
  double h_left;
  double h_right;
  double h_unperturbed;
  double h_perturbed;
  double u_left;
  double u_right;
  double bump_x_center;
  double bump_y_center;
  double bump_x_width;
  double bump_y_width;
  double bump_height;
  double perturbation_x_center;
  double perturbation_y_center;
  double perturbation_x_width;
  double perturbation_y_width;
};

#include "src/parameters/ShallowWaterProblemParameters.cc"

#endif
