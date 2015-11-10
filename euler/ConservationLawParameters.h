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
  /** number of components in solution */
  unsigned int n_components;

  // refinement parameters
  enum RefinementMode
  {
    space,
    time
  };
  RefinementMode refinement_mode;
  unsigned int initial_refinement_level;
  unsigned int n_refinement_cycles; /** number of refinement cycles */
  double refinement_fraction;
  double coarsening_fraction;

  enum TimeStepSizeMethod
  {
    constant_dt,
    cfl_condition
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

  enum TemporalDiscretization
  {
    SS,
    ERK1,
    ERK2,
    ERK3,
    ERK4,
    SDIRK22
  };
  TemporalDiscretization time_discretization;
  enum TemporalIntegrator
  {
    runge_kutta
  };
  TemporalIntegrator temporal_integrator;

  enum NonlinearSolverType
  {
    newton
  };
  NonlinearSolverType nonlinear_solver;

  enum LinearSolverType
  {
    direct,
    gmres,
    cg
  };
  LinearSolverType linear_solver;
  LinearSolverType mass_matrix_linear_solver;

  enum Verbosity
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

  static void declare_conservation_law_parameters(ParameterHandler & prm);
  void get_conservation_law_parameters(ParameterHandler & prm);

  unsigned int degree;
  unsigned int n_quadrature_points;

  enum ViscosityType
  {
    none,
    constant,
    old_first_order,
    max_principle,
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
};

#include "ConservationLawParameters.cc"

#endif
