#ifndef TransportParameters_cc
#define TransportParameters_cc

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/**
 * Parameters class for defining input parameters.
 */
template<int dim>
class TransportParameters
{
public:
  enum TemporalDiscretization
  {
    SS, FE, CN, BE, BDF2, SSP2, SSP3
  };

  enum RefinementMode
  {
    space, time
  };

  TransportParameters();
  static void declare_parameters(ParameterHandler &prm);
  void get_parameters(ParameterHandler &prm);

  unsigned int problem_id; // problem ID

  unsigned int degree; // polynomial degree of finite elements
  TemporalDiscretization time_discretization_option;
  RefinementMode refinement_mode;
  double time_refinement_factor; // reduction factor for time refinement
  bool use_adaptive_refinement; // option to use adaptive mesh refinement
  unsigned int n_refinement_cycles; // number of refinement cycles
  unsigned int initial_refinement_level; // initial level of refinement
  unsigned int linear_solver_option; // linear solver option
  unsigned int nonlinear_solver_option; // nonlinear solver option
  unsigned int nonlinear_tolerance; // nonlinear solver tolerance
  unsigned int viscosity_option; // option for viscosity used
  std::string entropy_string; // string for entropy function
  std::string entropy_derivative_string; // string for entropy derivative function
  double entropy_residual_coefficient; // value of entropy residual coefficient
  double jump_coefficient; // value of jump coefficient
  // temporal discretization for entropy residual
  TemporalDiscretization EV_time_discretization;
  bool output_mesh; // option to output mesh as .eps file
  bool is_steady_state; // is the problem steady-state?
  double end_time; // end time if transient problem is run
  double time_step_size; // time step size if transient problem is run
  bool output_exact_solution; // option to output exact solution
  unsigned int exact_solution_refinement_level; // refinement level exact solution
  bool output_initial_solution; // option to output initial solution
  bool output_DMP_bounds; // option to output DMP bounds
  double CFL_limit; // upper bound for the CFL number
  bool do_not_limit; // flag to turn off limiting for high order solution
  unsigned int n_quadrature_points; // number of quadrature points to use in formula
  bool save_convergence_results; // option to save convergence results
};

#include "TransportParameters.cc"
#endif
