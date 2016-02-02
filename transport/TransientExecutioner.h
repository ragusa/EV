#ifndef TransientExecutioner_cc
#define TransientExecutioner_cc

#include "Executioner.h"
#include "SSPRKTimeIntegrator.h"

using namespace dealii;

/**
 * Class for transient executioner.
 */
template <int dim>
class TransientExecutioner : public Executioner<dim>
{
public:
  /** \brief Alias for temporal discretization */
  using TemporalDiscretization =
    typename TransportParameters<dim>::TemporalDiscretization;

  /** \brief Alias for temporal discretization of entropy */
  using EntropyTemporalDiscretization =
    typename TransportParameters<dim>::EntropyTemporalDiscretization;

  TransientExecutioner(const TransportParameters<dim> & parameters,
                       Triangulation<dim> & triangulation,
                       const Tensor<1, dim> & transport_direction,
                       const FunctionParser<dim> & cross_section_function,
                       FunctionParser<dim> & source_function,
                       Function<dim> & incoming_function,
                       FunctionParser<dim> & initial_conditions_function,
                       const double & domain_volume,
                       PostProcessor<dim> & postprocessor,
                       const bool & source_is_time_dependent);

  void run() override;

private:
  void assembleMassMatrices();

  double enforceCFLCondition(const double & dt_proposed) const;

  void takeGalerkinStep_ssprk(SSPRKTimeIntegrator<dim> & ssprk);

  void takeGalerkinStep_theta(const double & dt, const double & t_new);

  void takeLowOrderStep_ssprk(SSPRKTimeIntegrator<dim> & ssprk);

  void takeLowOrderStep_theta(const double & dt, const double & t_new);

  void takeEntropyViscosityStep_ssprk(SSPRKTimeIntegrator<dim> & ssprk,
                                      EntropyViscosity<dim> & EV,
                                      const double & dt,
                                      const double & dt_old,
                                      const double & dt_older,
                                      const double & t_old);

  void takeEntropyViscosityFCTStep_ssprk(SSPRKTimeIntegrator<dim> & ssprk,
                                         FCT<dim> & fct,
                                         EntropyViscosity<dim> & EV,
                                         const double & dt,
                                         const double & dt_old);

  void takeGalerkinFCTStep_ssprk(SSPRKTimeIntegrator<dim> & ssprk,
                                 FCT<dim> & fct,
                                 const double & dt);

  const TemporalDiscretization temporal_discretization;

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

  Vector<double> ss_rhs_new;

  const double theta;

  Vector<double> tmp_vector;
};

#include "TransientExecutioner.cc"
#endif
