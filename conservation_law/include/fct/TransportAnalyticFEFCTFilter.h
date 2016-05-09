/**
 * \file TransportAnalyticFEFCTFilter.h
 * \brief Provides the header for the TransportAnalyticFEFCTFilter class.
 */
#ifndef TransportAnalyticFEFCTFilter_h
#define TransportAnalyticFEFCTFilter_h

#include "include/fct/ExplicitEulerFCTFilter.h"
#include "include/fct/TransportAnalyticSolutionBounds.h"

using namespace dealii;

/**
 * \brief Class for analytic transport bounds filter for forward Euler.
 */
template <int dim>
class TransportAnalyticFEFCTFilter : public FEFCTFilter<dim>
{
public:
  TransportAnalyticFEFCTFilter(
    const RunParameters & run_parameters,
    TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const QGauss<dim> & cell_quadrature,
    const std::shared_ptr<Limiter<dim>> limiter,
    const SparseMatrix<double> & lumped_mass_matrix,
    const std::map<unsigned int, double> & dirichlet_values);

  virtual void filter_antidiffusive_fluxes(
    const Vector<double> & old_solution,
    const double & dt,
    const Vector<double> & inviscid_ss_flux,
    const Vector<double> & ss_reaction,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const Vector<double> & ss_rhs,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix) override;

  bool check_bounds(const Vector<double> & new_solution) override;

  Vector<double> get_lower_solution_bound() const override;

  Vector<double> get_upper_solution_bound() const override;

protected:
  void compute_solution_bounds(const Vector<double> & old_solution,
                               const double & dt,
                               const Vector<double> & ss_reaction,
                               const Vector<double> & ss_rhs) override;

  /** \brief analytic transport solution bounds */
  TransportAnalyticSolutionBounds<dim> analytic_bounds;
};

#include "src/fct/TransportAnalyticFEFCTFilter.cc"

#endif
