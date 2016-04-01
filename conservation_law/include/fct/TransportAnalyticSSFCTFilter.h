/**
 * \file TransportAnalyticSSFCTFilter.h
 * \brief Provides the header for the TransportAnalyticSSFCTFilter class.
 */
#ifndef TransportAnalyticSSFCTFilter_h
#define TransportAnalyticSSFCTFilter_h

#include "include/fct/SteadyStateFCTFilter.h"
#include "include/fct/TransportAnalyticSolutionBounds.h"

using namespace dealii;

/**
 * \brief Class for analytic transport bounds filter for steady-state.
 */
template <int dim>
class TransportAnalyticSSFCTFilter : public SteadyStateFCTFilter<dim>
{
public:
  TransportAnalyticSSFCTFilter(
    const RunParameters & run_parameters,
    TransportProblemParameters<dim> & problem_parameters,
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

  bool check_bounds(const Vector<double> & new_solution) override;

  Vector<double> get_lower_solution_bound() const override;

  Vector<double> get_upper_solution_bound() const override;

protected:
  void compute_solution_bounds(const Vector<double> & solution,
                               const SparseMatrix<double> & low_order_ss_matrix,
                               const Vector<double> & ss_rhs) override;

  /** \brief analytic transport solution bounds */
  TransportAnalyticSolutionBounds<dim> analytic_bounds;
};

#include "src/fct/TransportAnalyticSSFCTFilter.cc"

#endif
