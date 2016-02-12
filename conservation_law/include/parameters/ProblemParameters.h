/**
 * \file ProblemParameters.h
 * \brief Provides the header for the ProblemParameters class.
 */
#ifndef ProblemParameters_h
#define ProblemParameters_h

#include <string>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/**
 * \class ProblemParameters
 * \brief Class for conservation law problem parameters
 */
template <int dim>
class ProblemParameters
{
public:
  ProblemParameters();

  static void declare_common_parameters(ParameterHandler & parameter_handler);

  void get_common_parameters(ParameterHandler & parameter_handler);

  bool valid_in_1d;
  bool valid_in_2d;
  bool valid_in_3d;

  std::string domain_shape;
  double x_start;
  double y_start;
  double z_start;
  double x_width;
  double y_width;
  double z_width;

  std::string boundary_conditions_type;
  bool use_exact_solution_as_dirichlet_bc;

  bool has_exact_solution;
  std::string exact_solution_type;

  bool has_default_end_time;
  double default_end_time;
};

#include "src/parameters/ProblemParameters.cc"

#endif
