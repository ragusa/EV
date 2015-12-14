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

  enum class ViscosityType
  {
    none,
    constant,
    low,
    DMP_low,
    DI_low,
    entropy
  };
  ViscosityType viscosity_type;
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
    algebraic,
    laplacian,
    graphtheoretic
  };
  DiffusionType diffusion_type;
};

#include "src/parameters/ConservationLawParameters.cc"

#endif
