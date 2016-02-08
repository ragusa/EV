#ifndef SteadyStateExecutioner_cc
#define SteadyStateExecutioner_cc

#include "Executioner.h"

using namespace dealii;

/**
 * \brief Class for steady-state executioner.
 */
template <int dim>
class SteadyStateExecutioner : public Executioner<dim>
{
public:
  /** \brief Alias for FCT initialization option */
  using FCTInitializationOption =
    typename TransportParameters<dim>::FCTInitializationOption;

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
  void compute_galerkin_solution();

  void compute_low_order_solution();

  void compute_entropy_viscosity_solution();

  void compute_FCT_solution();

  void compute_FCT_solution_cumulative_antidiffusion();
};

#include "SteadyStateExecutioner.cc"
#endif
