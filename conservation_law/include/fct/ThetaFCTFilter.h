/**
 * \file ThetaFCTFilter.h
 * \brief Provides the header for the ThetaFCTFilter class.
 */
#ifndef ThetaFCTFilter_h
#define ThetaFCTFilter_h

#include <deal.II/lac/sparse_matrix.h>

#include "include/fct/FCTFilter.h"

using namespace dealii;

/**
 * \brief Abstract base class for FCT filters for theta.
 */
template <int dim>
class ThetaFCTFilter : public FCTFilter<dim>
{
public:
  ThetaFCTFilter(const RunParameters & run_parameters,
                 const std::shared_ptr<Limiter<dim>> limiter,
                 const DoFHandler<dim> & dof_handler,
                 const FESystem<dim> & fe,
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
    SparseMatrix<double> & antidiffusion_matrix);

protected:
  virtual void compute_solution_bounds(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs_new,
    const Vector<double> & ss_rhs_old) = 0;

  virtual void compute_antidiffusion_bounds(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs_new,
    const Vector<double> & ss_rhs_old,
    const Vector<double> & cumulative_antidiffusion) = 0;

  /** \brief lumped mass matrix \f$\mathbf{M}^L\f$ */
  const SparseMatrix<double> * const lumped_mass_matrix;

  /** \brief theta parameter \f$\theta\f$ */
  const double theta;
};

#include "src/fct/ThetaFCTFilter.cc"

#endif
