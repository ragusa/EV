#ifndef TransportTransientExecutioner_cc
#define TransportTransientExecutioner_cc

#include "TransportExecutioner.h"
#include "SSPRKTimeIntegrator.h"
#include "RunParameters.h"

using namespace dealii;

/**
 * \brief Class for transient executioner of a transport problem.
 */
template <int dim>
class TransportTransientExecutioner : public TransportExecutioner<dim>
{
public:
  /** \brief Alias for temporal discretization classification */
  using TemporalDiscretizationClassification =
    typename RunParameters<dim>::TemporalDiscretizationClassification;

  TransportTransientExecutioner(
    const TransportRunParameters<dim> & parameters,
    Triangulation<dim> & triangulation,
    const Tensor<1, dim> & transport_direction,
    const double & transport_speed,
    const FunctionParser<dim> & cross_section_function,
    FunctionParser<dim> & source_function,
    Function<dim> & incoming_function,
    FunctionParser<dim> & initial_conditions_function,
    const double & domain_volume,
    PostProcessor<dim> & postprocessor,
    const bool & source_is_time_dependent,
    const double & nominal_dt);

  void run() override;

private:
  void assembleMassMatrices();

  double enforceCFLCondition(const double & dt_proposed) const;

  void compute_galerkin_solution_ssprk(SSPRKTimeIntegrator<dim> & ssprk);

  void compute_galerkin_solution_theta(const double & dt, const double & t_new);

  void compute_low_order_solution_ssprk(SSPRKTimeIntegrator<dim> & ssprk);

  void compute_low_order_solution_theta(const double & dt, const double & t_new);

  void compute_entropy_viscosity_solution_ssprk(SSPRKTimeIntegrator<dim> & ssprk,
                                                EntropyViscosity<dim> & EV,
                                                const double & dt_old,
                                                const double & dt_older,
                                                const double & t_old);

  void compute_entropy_viscosity_solution_theta(EntropyViscosity<dim> & EV,
                                                const double & dt,
                                                const double & dt_old,
                                                const double & t_new);

  void compute_entropy_viscosity_fct_solution_ssprk(
    SSPRKTimeIntegrator<dim> & ssprk,
    FCT<dim> & fct,
    EntropyViscosity<dim> & EV,
    const double & dt,
    const double & dt_old,
    const double & t_old);

  void compute_galerkin_fct_solution_ssprk(SSPRKTimeIntegrator<dim> & ssprk,
                                           FCT<dim> & fct,
                                           const double & dt,
                                           const double & t_old);

  void compute_fct_solution_theta(FCT<dim> & fct,
                                  const double & dt,
                                  const double & t_old);

  const TemporalDiscretizationClassification temporal_discretization;

  SparseMatrix<double> consistent_mass_matrix;

  SparseMatrix<double> lumped_mass_matrix;

  FunctionParser<dim> * const initial_conditions_function;

  double dt_nominal;

  double minimum_cell_diameter;

  Vector<double> old_solution;
  Vector<double> older_solution;
  Vector<double> oldest_solution;
  Vector<double> old_stage_solution;

  const bool source_is_time_dependent;

  SparseMatrix<double> high_order_diffusion_matrix_new;

  SparseMatrix<double> high_order_ss_matrix_new;

  Vector<double> ss_rhs_new;

  const double theta;

  Vector<double> tmp_vector;
};

#include "TransportTransientExecutioner.cc"
#endif
