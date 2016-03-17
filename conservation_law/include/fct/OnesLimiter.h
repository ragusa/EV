/**
 * \file OnesLimiter.h
 * \brief Provides the header for the OnesLimiter class.
 */
#ifndef OnesLimiter_h
#define OnesLimiter_h

#include "include/fct/Limiter.h"

using namespace dealii;

/**
 * \brief Class for limiter that produces no limitation.
 */
class OnesLimiter : public Limiter
{
public:
  OnesLimiter(const unsigned int & n_dofs);

  void compute_limiter_matrix(const SparseMatrix<double> & antidiffusion_matrix,
                              const DoFBounds & antidiffusion_bounds,
                              SparseMatrix<double> & limiter_matrix) override;
};

#include "src/fct/OnesLimiter.cc"

#endif
