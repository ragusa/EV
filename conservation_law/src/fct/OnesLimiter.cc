/**
 * \file OnesLimiter.cc
 * \brief Provides the function definitions for the OnesLimiter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 * \param[in] report_antidiffusion_  flag to report amount of accepted
 *            antidiffusion
 */
template <int dim>
OnesLimiter<dim>::OnesLimiter(const unsigned int & n_dofs_,
                              const bool & report_antidiffusion_)
  : Limiter<dim>(n_dofs_, report_antidiffusion_)
{
}

/**
 * \brief Computes the limiting coefficient matrix.
 *
 * \param[in] antidiffusion_matrix  matrix of antidiffusion fluxes
 *            \f$\mathbf{P}\f$
 * \param[in] antidiffusion_bounds  lower and upper bounds of antidiffusive
 *            fluxes into each node \f$\mathbf{Q}^-\f$ and \f$\mathbf{Q}^+\f$
 * \param[in] cumulative_antidiffusion_vector  antidiffusion accepted in previous
 *            iterations $\bar{\mathbf{p}}^{(\ell-1)}$
 * \param[out] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$
 */
template <int dim>
void OnesLimiter<dim>::compute_limiter_matrix(
  const SparseMatrix<double> &,
  const DoFBounds<dim> &,
  const Vector<double> &,
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
