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
template <int dim>
class OnesLimiter : public Limiter<dim>
{
public:
  OnesLimiter(const unsigned int & n_dofs,
              const bool & report_antidiffusion = false);

  void compute_limiter_matrix(
    const SparseMatrix<double> & antidiffusion_matrix,
    const DoFBounds<dim> & antidiffusion_bounds,
    const Vector<double> & cumulative_antidiffusion_vector,
    SparseMatrix<double> & limiter_matrix) override;
};

#include "src/fct/OnesLimiter.cc"

#endif
