/**
 * \file TransportAnalyticThetaFCTFilter.h
 * \brief Provides the header for the TransportAnalyticThetaFCTFilter class.
 */
#ifndef TransportAnalyticThetaFCTFilter_h
#define TransportAnalyticThetaFCTFilter_h

#include "include/fct/ThetaFCTFilter.h"
#include "include/fct/TransportAnalyticSolutionBounds.h"

using namespace dealii;

/**
 * \brief Class for analytic transport bounds filter for theta time
 *        discretizations.
 */
template <int dim>
class TransportAnalyticThetaFCTFilter : public ThetaFCTFilter<dim>
{
public:
  TransportAnalyticThetaFCTFilter(
    const RunParameters & run_parameters,
    TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const QGauss<dim> & cell_quadrature,
    const std::shared_ptr<Limiter<dim>> limiter,
                 const SparseMatrix<double> & lumped_mass_matrix,
                 const double & theta,
    const std::map<unsigned int, double> & dirichlet_values);

  virtual void filter_antidiffusive_fluxes(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs_new,
    const Vector<double> & ss_rhs_old,
    const Vector<double> & cumulative_antidiffusion,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix) override;

  bool check_bounds(const Vector<double> & new_solution) override;

  Vector<double> get_lower_solution_bound() const override;

  Vector<double> get_upper_solution_bound() const override;

protected:
  void compute_solution_bounds(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs_new,
    const Vector<double> & ss_rhs_old) override;

  /** \brief analytic transport solution bounds */
  TransportAnalyticSolutionBounds<dim> analytic_bounds;
};

#include "src/fct/TransportAnalyticThetaFCTFilter.cc"

#endif
