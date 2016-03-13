/**
 * \file ZeroesLimiter.h
 * \brief Provides the header for the ZeroesLimiter class.
 */
#ifndef ZeroesLimiter_h
#define ZeroesLimiter_h

#include <deal.II/lac/vector.h>

using namespace dealii;

/**
 * \brief Class for limiter that produces full limitation.
 */
class ZeroesLimiter : public Limiter
{
public:
  ZeroesLimiter(const unsigned int & n_dofs);

  void compute_limiter_matrix(const SparseMatrix<double> & antidiffusion_matrix,
                              const DoFBounds & antidiffusion_bounds,
                              SparseMatrix<double> & limiter_matrix) override;
};

#include "src/fct/ZeroesLimiter.cc"

#endif
