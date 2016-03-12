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
class Limiter
{
public:
  Limiter(const unsigned int & n_dofs);

  virtual void compute_limiter(const SparseMatrix<double> & antidiffusion_matrix,
    const DoFBounds & antidiffusion_bounds,
    SparseMatrix<double> & limiter_matrix) const = 0;

  void limit_antidiffusion(const SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix) const;

protected:
  /** \brief number of degrees of freedom */
  const unsigned int n_dofs;
};

#include "src/fct/Limiter.cc"

#endif
