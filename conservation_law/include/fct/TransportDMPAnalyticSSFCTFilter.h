/**
 * \file TransportDMPAnalyticSSFCTFilter.h
 * \brief Provides the header for the TransportDMPAnalyticSSFCTFilter class.
 */
#ifndef TransportDMPAnalyticSSFCTFilter_h
#define TransportDMPAnalyticSSFCTFilter_h

#include "include/fct/DMPSteadyStateFCTFilter.h"
#include "include/fct/TransportAnalyticSolutionBounds.h"

using namespace dealii;

/**
 * \brief Class for discrete maximum principle (DMP) with analytic transport
 *        bounds filter for steady-state.
 */
template <int dim>
class TransportDMPAnalyticSSFCTFilter : public DMPSteadyStateFCTFilter<dim>
{
public:
  TransportDMPAnalyticSSFCTFilter(
    const TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const QGauss<dim> & cell_quadrature,
    const std::shared_ptr<Limiter<dim>> limiter);

  virtual void filter_antidiffusive_fluxes(
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs,
    const Vector<double> & cumulative_antidiffusion,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix) override;

protected:
  /** \brief analytic transport solution bounds */
  TransportAnalyticSolutionBounds<dim> analytic_bounds;
};

#include "src/fct/TransportDMPAnalyticSSFCTFilter.cc"

#endif
