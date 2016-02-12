#ifndef TransportParameters_cc
#define TransportParameters_cc

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/**
 * Parameters class for defining input parameters.
 */
template <int dim>
class TransportParameters
{
public:
  TransportParameters();

  /** \brief Enumeration for types of temporal discretizations */
  enum class TemporalDiscretization
  {
    ss,
    theta,
    ssprk
  };
  TemporalDiscretization temporal_discretization;

  /** \brief Enumeration for types of temporal discretizations for entropy */
  enum class EntropyTemporalDiscretization
  {
    BE,
    CN,
    BDF2
  };
  EntropyTemporalDiscretization entropy_temporal_discretization;

  /** \brief Enumeration for SSPRK methods */
  enum class SSPRKMethod
  {
    FE,
    SSP2,
    SSP3
  };
  SSPRKMethod ssprk_method;

  /** \brief Enumeration for theta time discretization methods */
  enum class ThetaMethod
  {
    FE,
    CN,
    BE
  };
  ThetaMethod theta_method;

  /** \brief Theta parameter \f$\theta\f$ for theta time discretization */
  double theta;

  /** \brief Initialization option for FCT */
  enum class FCTInitializationOption
  {
    zero,
    low,
    high
  };
  FCTInitializationOption fct_initialization_option;

  static void declare_parameters(ParameterHandler & prm);

  void get_parameters(ParameterHandler & prm);

  // problem parameters
  unsigned int problem_id; // problem ID

  // time parameters
  double end_time;                 // end time if transient problem is run
  double time_step_size;           // time step size if transient problem is run
  bool use_adaptive_time_stepping; // option to adjust dt to satisfy CFL
  double CFL_limit;                // upper bound for the CFL number
  TemporalDiscretization time_discretization_option;

  // viscosity parameters
  unsigned int viscosity_option; // option for viscosity used

  // entropy viscosity parameters
  std::string entropy_string;            // string for entropy function
  std::string entropy_derivative_string; // string for entropy derivative function
  double entropy_residual_coefficient;   // value of entropy residual coefficient
  double jump_coefficient;               // value of jump coefficient

  // fct parameters
  bool do_not_limit; // flag to turn off limiting for high order solution
  bool skip_fct_if_bounds_satisfied; // option to skip FCT if high-order
                                     // solution already satisfies bounds
  bool use_cumulative_antidiffusion_algorithm; // option to use cumulative
                                               // antidiffusion algorithm for FCT
  bool include_analytic_bounds; // option to widen bounds with analytic bounds

  // refinement parameters
  bool refine_space;                     // option to refine space
  bool refine_time;                      // option to refine time
  double time_refinement_factor;         // reduction factor for time refinement
  bool use_adaptive_refinement;          // option to use adaptive mesh refinement
  unsigned int initial_refinement_level; // initial level of refinement
  unsigned int n_refinement_cycles;      // number of refinement cycles
  bool use_cell_size_for_convergence_rates; // use dx for rates, otherwise dt

  // finite element parameters
  unsigned int degree; // polynomial degree of finite elements
  unsigned int
    n_quadrature_points; // number of quadrature points to use in formula

  // linear solver options
  unsigned int linear_solver_option; // linear solver option

  // nonlinear solver options
  unsigned int nonlinear_solver_option; // nonlinear solver option
  double nonlinear_tolerance;           // nonlinear solver tolerance
  unsigned int nonlinear_max_iteration; // maximum iteration
  double relaxation_factor; // relaxation factor for nonlinear solution updates

  // output options
  bool output_mesh;           // option to output mesh as .eps file
  bool output_exact_solution; // option to output exact solution
  unsigned int exact_solution_refinement_level; // refinement level exact solution
  bool output_initial_solution;  // option to output initial solution
  bool output_DMP_bounds;        // option to output DMP bounds
  bool save_convergence_results; // option to save convergence results
  bool print_solution;           // option to print final solution
};

#include "TransportParameters.cc"
#endif
