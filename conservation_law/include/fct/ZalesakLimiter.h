/**
 * \file ZalesakLimiter.h
 * \brief Provides the header for the ZalesakLimiter class.
 */
#ifndef ZalesakLimiter_h
#define ZalesakLimiter_h

#include <deal.II/lac/vector.h>

#include "include/fct/Limiter.h"

using namespace dealii;

/**
 * \brief Class for the Zalesak FCT limiter.
 */
class ZalesakLimiter : public Limiter
{
public:
  ZalesakLimiter(const unsigned int & n_dofs);

  void compute_limiter_matrix(const SparseMatrix<double> & antidiffusion_matrix,
                              const DoFBounds & antidiffusion_bounds,
                              SparseMatrix<double> & limiter_matrix) override;

protected:
  /** \brief limiting coefficient for negative antidiffusive fluxes for a
             a node, without consideration of neighbor antidiffusion bounds,
             \f$\mathbf{L}^-\f$, commonly denoted as \f$\mathbf{R}^-\f$ */
  Vector<double> negative_limiter_vector;

  /** \brief limiting coefficient for positive antidiffusive fluxes for a
             a node, without consideration of neighbor antidiffusion bounds,
             \f$\mathbf{L}^+\f$, commonly denoted as \f$\mathbf{R}^+\f$ */
  Vector<double> positive_limiter_vector;
};

#include "src/fct/ZalesakLimiter.cc"

#endif
