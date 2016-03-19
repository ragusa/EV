#ifndef TransportSteadyStateExecutioner_cc
#define TransportSteadyStateExecutioner_cc

#include "include/executioners/TransportExecutioner.h"
#include "include/fct/TransportSteadyStateFCT.h"
#include "include/parameters/RunParameters.h"

using namespace dealii;

/**
 * \brief Class for steady-state executioner.
 */
template <int dim>
class TransportSteadyStateExecutioner : public TransportExecutioner<dim>
{
public:
  /** \brief Alias for scheme */
  using Scheme = typename TransportExecutioner<dim>::Scheme;
  /** \brief Alias for high-order scheme */
  using HighOrderScheme = typename TransportExecutioner<dim>::HighOrderScheme;
  /** \brief Alias for FCT initialization option */
  using FCTInitializationOption = typename RunParameters::FCTInitializationOption;

  TransportSteadyStateExecutioner(
    const TransportRunParameters & parameters,
    TransportProblemParameters<dim> & problem_parameters,
    Triangulation<dim> & triangulation,
    PostProcessor<dim> & postprocessor);

  void run() override;

protected:
  void compute_galerkin_solution();

  void compute_low_order_solution();

  void compute_high_order_solution();

  void compute_entropy_viscosity_solution();

  void compute_fct_solution();
};

#include "src/executioners/TransportSteadyStateExecutioner.cc"
#endif
