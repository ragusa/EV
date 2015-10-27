/**
 * \file ShallowWaterProblemParameters.h
 * \brief Provides the header for the ShallowWaterProblemParameters class.
 */
#ifndef ShallowWaterProblemParameters_h
#define ShallowWaterProblemParameters_h

#include <iostream>
#include <string>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/**
 * \class ShallowWaterProblemParameters
 * \brief Class for parameters related to the shallow water problems.
 */
template <int dim>
class ShallowWaterProblemParameters
{
public:
  ShallowWaterProblemParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  bool valid_in_1d;
  bool valid_in_2d;

  std::string domain_shape;
  double x_start;
  double y_start;
  double x_width;
  double y_width;

  std::string bathymetry_function;

  std::string boundary_conditions_type;
  bool use_exact_solution_as_dirichlet_bc;
  std::string dirichlet_function_height;
  std::string dirichlet_function_momentumx;
  std::string dirichlet_function_momentumy;

  std::string initial_conditions_height;
  std::string initial_conditions_momentumx;
  std::string initial_conditions_momentumy;

  bool has_exact_solution;
  std::string exact_solution_type;
  std::string exact_solution_height;
  std::string exact_solution_momentumx;
  std::string exact_solution_momentumy;

  bool has_default_end_time;
  double default_end_time;

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
  double perturbation_x_width;
};

#include "ShallowWaterProblemParameters.cc"

#endif
