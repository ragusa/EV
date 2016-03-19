/**
 * \file TransportSteadyStateFCT.h
 * \brief Provides the header for the TransportSteadyStateFCT class.
 */
#ifndef TransportSteadyStateFCT_h
#define TransportSteadyStateFCT_h

#include "include/fct/SteadyStateFCT.h"
#include "include/fct/TransportAnalyticFCTFilter.h"
#include "include/parameters/TransportRunParameters.h"
#include "include/parameters/TransportProblemParameters.h"

using namespace dealii;

/**
 * \brief Class for steady-state FCT for linear transport.
 */
template <int dim>
class TransportSteadyStateFCT : public SteadyStateFCT<dim>
{
public:
  TransportSteadyStateFCT(
    const TransportRunParameters & run_parameters,
    const TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler_);

protected:
  virtual std::shared_ptr<SteadyStateFCTFilter<dim>> create_filter(
    const std::string & filter_string) override;
};

#include "src/fct/TransportSteadyStateFCT.cc"

#endif
