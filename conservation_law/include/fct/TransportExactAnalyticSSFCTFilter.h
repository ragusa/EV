/**
 * \file TransportExactAnalyticSSFCTFilter.h
 * \brief Provides the header for the TransportExactAnalyticSSFCTFilter class.
 */
#ifndef TransportExactAnalyticSSFCTFilter_h
#define TransportExactAnalyticSSFCTFilter_h

#include "include/fct/TransportAnalyticSolutionBounds.h"
#include "include/fct/TransportExactAnalyticSSFCTFilter.h"

using namespace dealii;

/**
 * \brief Class to evaluate the analytic solution bounds \f$W_i^-\f$ and
 * \f$W_i^+\f$
 *        with the exact solution for steady-state.
 *
 * These solution bounds are not practical, as they require the exact solution,
 * but these bounds are useful as an academic exercise to determine if steady-
 * state FCT convergence difficulties are resolved when the solution bounds
 * are no longer implicit.
 */
template <int dim>
class TransportExactAnalyticSSFCTFilter : public TransportAnalyticSSFCTFilter<dim>
{
public:
  TransportExactAnalyticSSFCTFilter(
    const RunParameters & run_parameters,
    TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const QGauss<dim> & cell_quadrature,
    const std::shared_ptr<Limiter<dim>> limiter,
    const std::map<unsigned int, double> & dirichlet_values,
    const double & dx_min);

  void compute_solution_bounds(const Vector<double> & solution,
                               const SparseMatrix<double> & low_order_ss_matrix,
                               const Vector<double> & ss_rhs) override;

protected:
  /** \brief exact solution values at each node */
  Vector<double> exact_solution_values;
};

#include "src/fct/TransportExactAnalyticSSFCTFilter.cc"

#endif
