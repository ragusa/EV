/**
 * \file ZeroesLimiter.h
 * \brief Provides the header for the ZeroesLimiter class.
 */
#ifndef ZeroesLimiter_h
#define ZeroesLimiter_h

#include "include/fct/Limiter.h"

using namespace dealii;

/**
 * \brief Class for limiter that produces full limitation.
 */
template <int dim>
class ZeroesLimiter : public Limiter<dim>
{
public:
  ZeroesLimiter(const unsigned int & n_dofs);

  void compute_limiter_matrix(
    const SparseMatrix<double> & antidiffusion_matrix,
    const DoFBounds<dim> & antidiffusion_bounds,
    const Vector<double> & cumulative_antidiffusion_vector,
    SparseMatrix<double> & limiter_matrix) override;
};

#include "src/fct/ZeroesLimiter.cc"

#endif
