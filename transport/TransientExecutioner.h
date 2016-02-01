#ifndef TransientExecutioner_cc
#define TransientExecutioner_cc

#include "Executioner.h"
//#include "NonlinearSolver.h"
#include "SSPRKTimeIntegrator.h"

using namespace dealii;

/**
 * Class for transient executioner.
 */
template <int dim>
class TransientExecutioner : public Executioner<dim>
{
public:
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
  ~TransientExecutioner();

  void run() override;

private:
  void assembleMassMatrices();
  double enforceCFLCondition(const double & dt_proposed) const;

  void takeGalerkinStep(SSPRKTimeIntegrator<dim> & ssprk);
  void takeLowOrderStep(SSPRKTimeIntegrator<dim> & ssprk);
  void takeEntropyViscosityStep(SSPRKTimeIntegrator<dim> & ssprk,
                                EntropyViscosity<dim> & EV,
                                const double & dt,
                                const double & dt_old,
                                const double & dt_older,
                                const double & t_old);
  void takeEntropyViscosityFCTStep(SSPRKTimeIntegrator<dim> & ssprk,
                                   FCT<dim> & fct,
                                   EntropyViscosity<dim> & EV,
                                   const double & dt,
                                   const double & dt_old);
  void takeGalerkinFCTStep(SSPRKTimeIntegrator<dim> & ssprk,
                           FCT<dim> & fct,
                           const double & dt);

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
};

#include "TransientExecutioner.cc"
#endif
