/** \file ConservationLawParameters.h
 *  \brief Provides the header for the ConservationLawParameters class.
 */
#ifndef ConservationLawParameters_h
#define ConservationLawParameters_h

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/** \class ConservationLawParameters
 *  \brief Class for parameters related to solving conservation law equations.
 */
template <int dim>
class ConservationLawParameters
{
public:
  static void declare_conservation_law_parameters(ParameterHandler & prm);

  void get_conservation_law_parameters(ParameterHandler & prm);

  /** number of components in solution */
  unsigned int n_components;

  // refinement parameters
  enum class RefinementMode
  {
    space,
    time
  };
  RefinementMode refinement_mode;
  unsigned int initial_refinement_level;
  unsigned int n_refinement_cycles; /** number of refinement cycles */
  double refinement_fraction;
  double coarsening_fraction;

  enum class TimeStepSizeMethod
  {
    constant,
    cfl
  };
  TimeStepSizeMethod time_step_size_method;
  bool use_default_end_time;
  double end_time;
  double time_step_size;
  double cfl;
  double steady_state_tolerance;

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

  enum class TemporalDiscretization
  {
    FE,
    SSP2,
    SSP3
  };
  TemporalDiscretization time_discretization;
  enum class TemporalIntegrator
  {
    runge_kutta
  };
  TemporalIntegrator temporal_integrator;

  enum class NonlinearSolverType
  {
    newton
  };
  NonlinearSolverType nonlinear_solver;

  enum class LinearSolverType
  {
    direct,
    gmres,
    cg
  };
  LinearSolverType linear_solver;
  LinearSolverType mass_matrix_linear_solver;

  enum class Verbosity
  {
    quiet,
    verbose
  };
  Verbosity nonlinear_verbosity;
  Verbosity linear_verbosity;

  double linear_atol;
  double linear_rtol;
  unsigned int max_linear_iterations;

  double nonlinear_atol;
  double nonlinear_rtol;
  unsigned int max_nonlinear_iterations;
  double damping;

  unsigned int degree;
  unsigned int n_quadrature_points;

  enum class Scheme
  {
    low,
    high,
    fct
  };
  Scheme scheme;

  enum class LowOrderScheme
  {
    constant,
    standard,
    dmp,
    di_visc,
    di_diff
  };
  LowOrderScheme low_order_scheme;

  enum class HighOrderScheme
  {
    galerkin,
    entropy_visc,
    entropy_diff
  };
  HighOrderScheme high_order_scheme;

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
  std::string entropy_viscosity_smoothing;
  unsigned int entropy_viscosity_smoothing_weight;

  enum class DiffusionType
  {
    none,
    laplacian,
    graphtheoretic,
    DI,
    entropy
  };

  enum class AntidiffusionType
  {
    limited,
    full,
    none
  };
  AntidiffusionType antidiffusion_type;

  enum class FCTSynchronizationType
  {
    none,
    min,
    compound
  };
  FCTSynchronizationType fct_synchronization_type;

  bool use_star_states_in_fct_bounds;
};

#include "src/parameters/ConservationLawParameters.cc"

#endif
