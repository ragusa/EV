/**
 * \file TransportUpwindAnalyticFEFCTFilter.h
 * \brief Provides the header for the TransportUpwindAnalyticFEFCTFilter class.
 */
#ifndef TransportUpwindAnalyticFEFCTFilter_h
#define TransportUpwindAnalyticFEFCTFilter_h

#include "include/fct/ExplicitEulerFCTFilter.h"
#include "include/fct/TransportUpwindAnalyticSolutionBounds.h"

using namespace dealii;

/**
 * \brief Class for upwind analytic transport bounds filter for forward Euler.
 */
template <int dim>
class TransportUpwindAnalyticFEFCTFilter : public ExplicitEulerFCTFilter<dim>
{
public:
  TransportUpwindAnalyticFEFCTFilter(
    const RunParameters & run_parameters,
    TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
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
    const double & t_old,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix) override;

  bool check_bounds(const Vector<double> & new_solution) override;

  Vector<double> get_lower_solution_bound() const override;

  Vector<double> get_upper_solution_bound() const override;

protected:
  void compute_solution_bounds(const Vector<double> & old_solution,
                               const double & dt,
                               const Vector<double> & ss_reaction,
                               const Vector<double> & ss_rhs,
                               const double & t_old) override;

  /** \brief analytic transport solution bounds */
  TransportUpwindAnalyticSolutionBounds<dim> analytic_bounds;

  /** \brief transport speed */
  const double speed;
};

#include "src/fct/TransportUpwindAnalyticFEFCTFilter.cc"

#endif
