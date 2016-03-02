#ifndef TransportSteadyStateExecutioner_cc
#define TransportSteadyStateExecutioner_cc

#include "TransportExecutioner.h"
#include "RunParameters.h"

using namespace dealii;

/**
 * \brief Class for steady-state executioner.
 */
template <int dim>
class TransportSteadyStateExecutioner : public TransportExecutioner<dim>
{
public:
  /** \brief Alias for FCT initialization option */
  using FCTInitializationOption =
    typename RunParameters<dim>::FCTInitializationOption;

  TransportSteadyStateExecutioner(
    const TransportRunParameters<dim> & parameters,
    Triangulation<dim> & triangulation,
    const Tensor<1, dim> & transport_direction,
    const double & transport_speed,
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
};

#include "TransportSteadyStateExecutioner.cc"
#endif
