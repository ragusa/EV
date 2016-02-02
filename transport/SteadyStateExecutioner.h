#ifndef SteadyStateExecutioner_cc
#define SteadyStateExecutioner_cc

#include "Executioner.h"

using namespace dealii;

/**
 * Class for steady-state executioner.
 */
template <int dim>
class SteadyStateExecutioner : public Executioner<dim>
{
public:
  SteadyStateExecutioner(const TransportParameters<dim> & parameters,
                         Triangulation<dim> & triangulation,
                         const Tensor<1, dim> & transport_direction,
                         const FunctionParser<dim> & cross_section_function,
                         FunctionParser<dim> & source_function,
                         Function<dim> & incoming_function,
                         const double & domain_volume,
                         PostProcessor<dim> & postprocessor);

  void run() override;

protected:
  void compute_FCT_solution();
};

#include "SteadyStateExecutioner.cc"
#endif
