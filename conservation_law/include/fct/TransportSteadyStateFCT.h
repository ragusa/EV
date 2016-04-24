/**
 * \file TransportSteadyStateFCT.h
 * \brief Provides the header for the TransportSteadyStateFCT class.
 */
#ifndef TransportSteadyStateFCT_h
#define TransportSteadyStateFCT_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include "include/fct/SteadyStateFCT.h"
#include "include/fct/TransportAnalyticSSFCTFilter.h"
#include "include/fct/TransportDMPAnalyticSSFCTFilter.h"
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
  TransportSteadyStateFCT(const TransportRunParameters & run_parameters,
                          TransportProblemParameters<dim> & problem_parameters,
                          const DoFHandler<dim> & dof_handler,
                          const FESystem<dim> & fe,
                          const std::map<unsigned int, double> & dirichlet_values,
                          const QGauss<dim> & cell_quadrature);

protected:
  virtual std::shared_ptr<SteadyStateFCTFilter<dim>> create_filter(
    const std::string & filter_string) override;

  /** \brief problem parameters */
  TransportProblemParameters<dim> * const problem_parameters;

  /** \brief cell quadrature */
  const QGauss<dim> * const cell_quadrature;
};

#include "src/fct/TransportSteadyStateFCT.cc"

#endif
