/**
 * \file OnesLimiter.cc
 * \brief Provides the function definitions for the OnesLimiter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 */
OnesLimiter::OnesLimiter(const unsigned int & n_dofs_) : Limiter(n_dofs_) {}
/**
 * \brief Computes the limiting coefficient matrix.
 *
 * \param[in] antidiffusion_matrix  matrix of antidiffusion fluxes
 *            \f$\mathbf{P}\f$
 * \param[in] antidiffusion_bounds  lower and upper bounds of antidiffusive
 *            fluxes into each node \f$\mathbf{Q}^-\f$ and \f$\mathbf{Q}^+\f$
 * \param[out] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$
 */
void OnesLimiter::compute_limiter_matrix(const SparseMatrix<double> &,
                                         const DoFBounds &,
                                         SparseMatrix<double> & limiter_matrix)
{
  // set all entries to one
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    SparseMatrix<double>::iterator it = limiter_matrix.begin(i);
    SparseMatrix<double>::iterator it_end = limiter_matrix.end(i);
    for (; it != it_end; ++it)
      it->value() = 1.0;
  }
}
