/**
 * \file ZeroesLimiter.cc
 * \brief Provides the function definitions for the ZeroesLimiter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 */
ZeroesLimiter::ZeroesLimiter(const unsigned int & n_dofs_) : Limiter(n_dofs_) {}
/**
 * \brief Computes the limiting coefficient matrix.
 *
 * \param[in] antidiffusion_matrix  matrix of antidiffusion fluxes
 *            \f$\mathbf{P}\f$
 * \param[in] antidiffusion_bounds  lower and upper bounds of antidiffusive
 *            fluxes into each node \f$\mathbf{Q}^-\f$ and \f$\mathbf{Q}^+\f$
 * \param[out] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$
 */
void ZeroesLimiter::compute_limiter_matrix(const SparseMatrix<double> &,
                                           const DoFBounds &,
                                           SparseMatrix<double> & limiter_matrix)
{
  // set all entries to zero
  limiter_matrix = 0;
}
