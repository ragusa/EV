/**
 * \file Limiter.h
 * \brief Provides the header for the Limiter class.
 */
#ifndef Limiter_h
#define Limiter_h

#include <deal.II/lac/sparse_matrix.h>

#include "include/fct/DoFBounds.h"

using namespace dealii;

/**
 * \brief Abstract base class for FCT limiters.
 */
template <int dim>
class Limiter
{
public:
  Limiter(const unsigned int & n_dofs);

  void compute_limiter_matrix(const SparseMatrix<double> & antidiffusion_matrix,
                              const DoFBounds<dim> & antidiffusion_bounds,
                              SparseMatrix<double> & limiter_matrix);

  virtual void compute_limiter_matrix(
    const SparseMatrix<double> & antidiffusion_matrix,
    const DoFBounds<dim> & antidiffusion_bounds,
    const Vector<double> & cumulative_antidiffusion_vector,
    SparseMatrix<double> & limiter_matrix) = 0;

  void apply_limiter_matrix(const SparseMatrix<double> & limiter_matrix,
                            SparseMatrix<double> & antidiffusion_matrix) const;

protected:
  /** \brief number of degrees of freedom */
  const unsigned int n_dofs;

  /** \brief zero vector, to be used as a default for cumulative antidiffusion */
  const Vector<double> zero_vector;
};

#include "src/fct/Limiter.cc"

#endif
