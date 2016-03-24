/**
 * \file ExplicitEulerFCTFilter.h
 * \brief Provides the header for the ExplicitEulerFCTFilter class.
 */
#ifndef ExplicitEulerFCTFilter_h
#define ExplicitEulerFCTFilter_h

#include <deal.II/lac/sparse_matrix.h>

#include "include/fct/FCTFilter.h"

using namespace dealii;

/**
 * \brief Abstract base class for FCT filters for explicit Euler.
 */
template <int dim>
class ExplicitEulerFCTFilter : public FCTFilter<dim>
{
public:
  ExplicitEulerFCTFilter(const RunParameters & run_parameters,
                         const std::shared_ptr<Limiter<dim>> limiter,
                         const DoFHandler<dim> & dof_handler,
                         const FESystem<dim> & fe,
                         const SparseMatrix<double> & lumped_mass_matrix);

  virtual void filter_antidiffusive_fluxes(
    const Vector<double> & old_solution,
    const double & dt,
    const Vector<double> & inviscid_ss_flux,
    const Vector<double> & ss_reaction,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const Vector<double> & ss_rhs,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix);

protected:
  virtual void compute_solution_bounds(const Vector<double> & old_solution,
                                       const double & dt,
                                       const Vector<double> & ss_reaction,
                                       const Vector<double> & ss_rhs) = 0;

  virtual void compute_antidiffusion_bounds(
    const Vector<double> & old_solution,
    const double & dt,
    const Vector<double> & inviscid_ss_flux,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const Vector<double> & ss_rhs) = 0;

  /** \brief lumped mass matrix \f$\mathbf{M}^L\f$ */
  const SparseMatrix<double> * const lumped_mass_matrix;
};

#include "src/fct/ExplicitEulerFCTFilter.cc"

#endif
