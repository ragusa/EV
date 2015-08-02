#ifndef TransientExecutioner_cc
#define TransientExecutioner_cc

#include "Executioner.h"
#include "SSPRKTimeIntegrator.h"

using namespace dealii;

/**
 * Class for transient executioner.
 */
template<int dim>
class TransientExecutioner : public Executioner<dim>
{
public:

  TransientExecutioner(const TransportParameters<dim> & parameters,
    const Triangulation<dim> & triangulation,
    const Tensor<1, dim> & transport_direction,
    const FunctionParser<dim> & cross_section_function,
    FunctionParser<dim> & source_function, Function<dim> & incoming_function,
    FunctionParser<dim> & initial_conditions_function,
    PostProcessor<dim> & postprocessor,
    const bool & source_is_time_dependent);
  ~TransientExecutioner();

  void run() override;

private:

  void assembleMassMatrices();
  double enforceCFLCondition(double & dt);

  void takeGalerkinStep(SSPRKTimeIntegrator<dim> & ssprk);
  void takeLowOrderStep(SSPRKTimeIntegrator<dim> & ssprk);

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
