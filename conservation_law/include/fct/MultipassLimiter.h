/**
 * \file MultipassLimiter.h
 * \brief Provides the header for the MultipassLimiter class.
 */
#ifndef MultipassLimiter_h
#define MultipassLimiter_h

#include "include/fct/Limiter.h"

using namespace dealii;

/**
 * \brief Class for multi-pass limiter.
 *
 * This limiter allows multiple sweeps through any other limiter to allow
 * more antidiffusion to be kept.
 */
template <int dim>
class MultipassLimiter : public Limiter<dim>
{
public:
  MultipassLimiter(const unsigned int & n_dofs,
                   std::shared_ptr < Limiter<dim> limiter);

  void compute_limiter_matrix(
    const SparseMatrix<double> & antidiffusion_matrix,
    const DoFBounds<dim> & antidiffusion_bounds,
    const Vector<double> & cumulative_antidiffusion_vector,
    SparseMatrix<double> & limiter_matrix) override;

  void apply_limiter_matrix(const SparseMatrix<double> & limiter_matrix,
                            SparseMatrix<double> & antidiffusion_matrix) const;

protected:
  /** \brief number of degrees of freedom */
  const unsigned int n_dofs;

  /** \brief limiter for which multiple passes are to be made */
  std::shared_ptr<Limiter<dim>> limiter;
};

#include "src/fct/MultipassLimiter.cc"

#endif
