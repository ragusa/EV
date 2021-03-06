/**
 * \file SteadyStateFCTFilter.h
 * \brief Provides the header for the SteadyStateFCTFilter class.
 */
#ifndef SteadyStateFCTFilter_h
#define SteadyStateFCTFilter_h

#include <deal.II/lac/sparse_matrix.h>

#include "include/fct/FCTFilter.h"

using namespace dealii;

/**
 * \brief Abstract base class for FCT filters for steady-state.
 */
template <int dim>
class SteadyStateFCTFilter : public FCTFilter<dim>
{
public:
  SteadyStateFCTFilter(const RunParameters & run_parameters,
                       const std::shared_ptr<Limiter<dim>> limiter,
                       const DoFHandler<dim> & dof_handler,
                       const FESystem<dim> & fe,
                       const std::map<unsigned int, double> & dirichlet_values);

  virtual void filter_antidiffusive_fluxes(
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs,
    const Vector<double> & cumulative_antidiffusion,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix);

  virtual void compute_solution_bounds(
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs) = 0;

protected:
  void compute_antidiffusion_bounds(
    const DoFBounds<dim> & solution_bounds,
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs,
    const Vector<double> & cumulative_antidiffusion);
};

#include "src/fct/SteadyStateFCTFilter.cc"

#endif
