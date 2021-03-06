/** \file RunParameters.h
 *  \brief Provides the header for the RunParameters class.
 */
#ifndef RunParameters_h
#define RunParameters_h

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/**
 * \brief Class for run parameters.
 */
// template <int dim>
class RunParameters
{
public:
  static void declare_run_parameters(ParameterHandler & prm);

  void get_run_parameters(ParameterHandler & prm);

  /* -----------------------------------------------------------------------
   *                              problem name
   * ----------------------------------------------------------------------- */
  /** \brief name of the problem to be run */
  std::string problem_name;

  /* -----------------------------------------------------------------------
   *                              scheme
   * ----------------------------------------------------------------------- */
  /** \brief Enumeration for which scheme to run */
  enum class Scheme
  {
    low,
    high,
    fct
  };
  /** \brief Enumeration for which low-order scheme to run */
  enum class LowOrderScheme
  {
    constant,
    lax,
    dmp,
    di_visc,
    di_diff
  };
  /** \brief Enumeration for which high-order scheme to run */
  enum class HighOrderScheme
  {
    galerkin,
    entropy_visc,
    entropy_diff
  };

  /** \brief scheme to run */
  Scheme scheme;
  /** \brief low-order scheme to run */
  LowOrderScheme low_order_scheme;
  /** \brief high-order scheme to run */
  HighOrderScheme high_order_scheme;

  /* -----------------------------------------------------------------------
   *                              refinement
   * ----------------------------------------------------------------------- */
  /** \brief number of refinement cycles */
  unsigned int n_refinement_cycles;
  /** \brief option to refine space */
  bool refine_space;
  /** \brief option to refine time */
  bool refine_time;
  /** \brief initial refinement level */
  unsigned int initial_refinement_level;
  /** \brief option to use adaptive mesh refinement */
  bool use_adaptive_mesh_refinement;
  /** \brief reduction factor if refining time */
  double time_refinement_factor;
  /** \brief fraction of cells to refine if using adaptive mesh refinement */
  double refinement_fraction;
  /** \brief fraction of cells to coarsen if using adaptive mesh refinement */
  double coarsening_fraction;
  /** \brief flag to use dx for rates, otherwise dt */
  bool use_cell_size_for_convergence_rates;

  /* -----------------------------------------------------------------------
   *                              time
   * ----------------------------------------------------------------------- */
  /** \brief Enumeration for how to determine time step size */
  enum class TimeStepSizeOption
  {
    constant,
    cfl,
    cfl_dmp,
    cfl_di
  };

  /** \brief option for how to determine time step size */
  TimeStepSizeOption time_step_size_option;
  /** \brief CFL limit to impose if chose to compute time step size from CFL */
  double cfl;
  /** \brief time step size to use if chose constant time step size */
  double time_step_size;
  /** \brief option to use default end time for problem if available */
  bool use_default_end_time;
  /** \brief end time if default end time is not used */
  double end_time;
  /** \brief steady-state tolerance */
  double steady_state_tolerance;

  /* -----------------------------------------------------------------------
   *                              temporal discretization
   * ----------------------------------------------------------------------- */
  /** \brief Enumeration for classification of temporal discretizations */
  enum class TemporalDiscretizationClassification
  {
    ss,
    theta,
    ssprk
  };
  /** \brief Enumeration for SSPRK discretizations */
  enum class SSPRKDiscretization
  {
    FE,
    SSP2,
    SSP3
  };
  /** \brief Enumeration for theta discretizations */
  enum class ThetaDiscretization
  {
    FE,
    CN,
    BE
  };

  /** \brief classification of temporal discretization to use */
  TemporalDiscretizationClassification temporal_discretization;
  /** \brief SSPRK discretization to use */
  SSPRKDiscretization ssprk_discretization;
  /** \brief theta discretization to use */
  ThetaDiscretization theta_discretization;
  /** \brief Theta parameter \f$\theta\f$ for theta time discretization */
  double theta;

  /* -----------------------------------------------------------------------
   *                              artificial viscosity
   * ----------------------------------------------------------------------- */
  /** \brief Enumeration for type of artificial viscosity */
  enum class ViscosityType
  {
    none,
    constant,
    lax,
    DMP,
    DI,
    entropy,
    high
  };
  /** \brief Enumeration for type of artificial diffusion */
  enum class DiffusionType
  {
    none,
    laplacian,
    graphtheoretic,
    DI,
    entropy
  };
  /** \brief Enumeration for type of smoothing for entropy viscosity */
  enum class EntropyViscositySmoothingOption
  {
    none,
    max,
    average
  };

  /** \brief Coefficient for entropy residual in entropy viscosity */
  double entropy_residual_coef;
  /** \brief Coefficient for entropy jump in entropy viscosity */
  double entropy_jump_coef;
  /** \brief smoothing option for entropy viscosity */
  EntropyViscositySmoothingOption entropy_viscosity_smoothing_option;
  /** \brief smoothing weight for entropy viscosity */
  unsigned int entropy_viscosity_smoothing_weight;
  /** \brief value for constant viscosity */
  double constant_viscosity_value;
  /** \brief Coefficient for Lax viscosity */
  double lax_viscosity_coef;
  /** \brief flag to use low-order viscosity for first time step */
  bool use_low_order_viscosity_for_first_time_step;

  /* -----------------------------------------------------------------------
   *                              fct
   * ----------------------------------------------------------------------- */
  /** \brief Enumeration for limiter options */
  enum class LimiterOption
  {
    ones,
    zeroes,
    zalesak
  };
  /** \brief Enumeration for type of synchronization to apply for FCT */
  enum class FCTSynchronizationType
  {
    none,
    min,
    compound
  };
  /** \brief Enumeration for initialization option for implicit/SS FCT */
  enum class FCTInitializationOption
  {
    zero,
    low,
    high
  };

  /** \brief string of sequence of FCT filters */
  std::string filter_sequence_string;
  /** \brief limiter option */
  LimiterOption limiter_option;
  /** \brief flag to use multi-pass limiting */
  bool use_multipass_limiting;
  /** \brief percent tolerance for multi-pass limiting */
  double multipass_limiting_percent_tolerance;
  /** \brief max value for Dirichlet node limiting coefficients */
  double dirichlet_limiting_coefficient;
  /** \brief option to force correct signs of antidiffusion bounds */
  bool enforce_antidiffusion_bounds_signs;
  /** \brief type of synchronization to apply for FCT */
  FCTSynchronizationType fct_synchronization_type;
  /** \brief Initialization option for implicit and steady-state FCT */
  FCTInitializationOption fct_initialization_option;
  /** \brief option to skip FCT if high-order solution satisfies bounds */
  bool skip_fct_if_bounds_satisfied;
  /** \brief option to use the cumulative antidiffusion FCT algorithm */
  bool use_cumulative_antidiffusion_algorithm;
  /** \brief Flag to include star states in FCT bounds */
  bool use_star_states_in_fct_bounds;
  /** \brief number of sampling points for min/max in upwind solution bounds */
  unsigned int upwind_bounds_sampling_points;
  /** \brief Flag to check FCT bounds */
  bool check_fct_bounds;
  /** \brief Flag to output the limiter matrix */
  bool output_limiter_matrix;
  /** \brief Flag to output the final FCT bounds */
  bool output_final_fct_bounds;
  /** \brief Flag to output the transient FCT bounds */
  bool output_transient_fct_bounds;

  /* -----------------------------------------------------------------------
   *                              linear solver
   * ----------------------------------------------------------------------- */
  /** \brief Enumeration for types of linear solvers */
  enum class LinearSolverType
  {
    direct,
    gmres
  };
  /// preconditioner types
  enum class PreconditionerType
  {
    none,
    jacobi,
    sor,
    ssor
  };

  /** \brief type of linear solver to use */
  LinearSolverType linear_solver_type;
  /// preconditioner type
  PreconditionerType preconditioner_type;
  /** \brief maximum number of linear iterations */
  unsigned int max_linear_iterations;
  /** \brief linear solve tolerance */
  double linear_tolerance;
  /// relaxation parameter for preconditioner
  double preconditioner_relaxation;
  /// print linear residuals
  bool print_linear_residuals;

  /* -----------------------------------------------------------------------
   *                              nonlinear solver
   * ----------------------------------------------------------------------- */
  /** \brief nonlinear tolerance */
  double nonlinear_tolerance;
  /** \brief maximum number of nonlinear iterations */
  unsigned int nonlinear_max_iterations;
  /** \brief relaxation factor for damping solution updates */
  double relaxation_factor;

  /* -----------------------------------------------------------------------
   *                              finite element
   * ----------------------------------------------------------------------- */
  /** \brief finite element degree */
  unsigned int degree;

  /* -----------------------------------------------------------------------
   *                              quadrature
   * ----------------------------------------------------------------------- */
  /** \brief number of quadrature points */
  unsigned int n_quadrature_points;

  /* -----------------------------------------------------------------------
   *                              output
   * ----------------------------------------------------------------------- */
  /** \brief Verbosity level (0 is least verbose) */
  unsigned int verbosity_level;
  /** \brief Period of transient output (0 = no transient, 1 = every step) */
  unsigned int output_period;
  /** \brief Limit of size for transient output */
  double transient_output_size_limit;
  /** \brief Flag to output mesh */
  bool output_mesh;
  /** \brief Flag to output mass matrix */
  bool output_mass_matrix;
  /** \brief Flag to output viscosity */
  bool output_viscosity;
  /** \brief Flag to output viscosity transient */
  bool output_viscosity_transient;
  /** \brief Flag to output exact solution */
  bool output_exact_solution;
  /** \brief refinement level of exact solution output */
  unsigned int exact_solution_refinement_level;
  /** \brief Flag to save convergence results to file */
  bool save_convergence_results;
  /** \brief Flag to print final solution to screen */
  bool print_final_solution;
  /** \brief Output directory (relative to executable directory) */
  std::string output_directory;
  /** \brief Flag to use a subdirectory named by problem name */
  bool use_problem_name_output_subdirectory;
  /** \brief Flag to append the scheme to the output filename */
  bool append_scheme_to_output_filename;
  /** \brief Flag to append the time discretization to the output filename */
  bool append_time_discretization_to_output_filename;
};

#include "src/parameters/RunParameters.cc"

#endif
