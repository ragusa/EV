/**
 * \file TransportUpwindAnalyticSSFCTFilter.h
 * \brief Provides the header for the TransportUpwindAnalyticSSFCTFilter class.
 */
#ifndef TransportUpwindAnalyticSSFCTFilter_h
#define TransportUpwindAnalyticSSFCTFilter_h

#include "include/fct/SteadyStateFCTFilter.h"
#include "include/fct/TransportUpwindAnalyticSolutionBounds.h"

using namespace dealii;

/**
 * \brief Class for analytic transport bounds filter for steady-state.
 */
template <int dim>
class TransportUpwindAnalyticSSFCTFilter : public SteadyStateFCTFilter<dim>
{
public:
  TransportUpwindAnalyticSSFCTFilter(
    const RunParameters & run_parameters,
    TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const std::shared_ptr<Limiter<dim>> limiter,
    const std::map<unsigned int, double> & dirichlet_values,
    const double & dx_min);

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

  virtual void compute_solution_bounds(
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs) override;

protected:
  /** \brief analytic transport solution bounds */
  TransportUpwindAnalyticSolutionBounds<dim> analytic_bounds;

  /** \brief minimum cell diameter */
  const double dx_min;
};

#include "src/fct/TransportUpwindAnalyticSSFCTFilter.cc"

#endif
