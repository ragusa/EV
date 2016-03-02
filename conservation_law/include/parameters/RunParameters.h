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
template <int dim>
class RunParameters
{
public:
  static void declare_run_parameters(ParameterHandler & prm);

  void get_run_parameters(ParameterHandler & prm);

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
  bool use_adaptive_refinement;
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
    cfl
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
  TemporalDiscretization temporal_discretization;
  /** \brief SSPRK discretization to use */
  SSPRKDiscretization ssprk_discretization;
  /** \brief theta discretization to use */
  ThetaMethod theta_discretization;
  /** \brief Theta parameter \f$\theta\f$ for theta time discretization */
  double theta;

  /* -----------------------------------------------------------------------
   *                              linear solver
   * ----------------------------------------------------------------------- */
  enum class LinearSolverType
  {
    direct
  };
  LinearSolverType linear_solver;
  LinearSolverType mass_matrix_linear_solver;

  /* -----------------------------------------------------------------------
   *                              nonlinear solver
   * ----------------------------------------------------------------------- */
  unsigned int max_nonlinear_iterations;
  double damping;



  enum class ViscosityType
  {
    none,
    constant,
    low,
    DMP,
    DI,
    entropy,
    high
  };

  double constant_viscosity_value;
  double first_order_viscosity_coef;
  bool use_low_order_viscosity_for_first_time_step;
  double entropy_residual_coef;
  double entropy_jump_coef;
  enum class EntropyViscositySmoothingOption
  {
    none,
    max,
    average
  };
  EntropyViscositySmoothingOption entropy_viscosity_smoothing_option;
  unsigned int entropy_viscosity_smoothing_weight;

  enum class DiffusionType
  {
    none,
    laplacian,
    graphtheoretic,
    DI,
    entropy
  };

  /** \brief Enumeration for type of bounds to impose on FCT solution */
  enum class FCTBoundsType
  {
    dmp
  };
  /** \brief Type of bounds to impose on FCT solution */
  FCTBoundsType fct_bounds_type;

  enum class AntidiffusionOption
  {
    limited,
    full,
    none
  };
  AntidiffusionOption antidiffusion_option;

  enum class FCTSynchronizationType
  {
    none,
    min,
    compound
  };
  FCTSynchronizationType fct_synchronization_type;

  /** \brief Enumeration for set of variables to limit in FCT */
  enum class FCTLimitationType
  {
    conservative,
    characteristic,
    sw_heightonly
  };
  /** \brief Set of variables to limit in FCT */
  FCTLimitationType fct_limitation_type;

  /** \brief Initialization option for implicit and steady-state FCT */
  enum class FCTInitializationOption
  {
    zero,
    low,
    high
  };
  FCTInitializationOption fct_initialization_option;

  bool use_star_states_in_fct_bounds;

  bool output_limiter_matrix;

  bool output_fct_bounds;

  unsigned int output_period;
  double max_transient_output_size;
  bool output_mesh;
  bool output_mass_matrix;
  bool output_viscosity;
  bool output_viscosity_transient;
  bool output_exact_solution;
  unsigned int exact_solution_refinement_level;
  bool save_convergence_results;
  bool print_final_solution;
  unsigned int verbosity_level;
};

#include "src/parameters/RunParameters.cc"

#endif
