/**
 * \file MultipassLimiter.cc
 * \brief Provides the function definitions for the MultipassLimiter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 * \param[in] limiter_  limiter for which multiple passes are to be made
 */
template <int dim>
MultipassLimiter<dim>::MultipassLimiter(const unsigned int & n_dofs_,
                                        std::shared_ptr < Limiter<dim> limiter_)
  : Limiter<dim>(n_dofs_), limiter(limiter_)
{
}

/**
 * \brief Computes the limiting coefficient matrix.
 *
 * This function takes passes through a limiter until no more antidiffusion
 * is accepted.
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
void MultipassLimiter<dim>::compute_limiter_matrix(
  const SparseMatrix<double> & antidiffusion_matrix,
  const DoFBounds<dim> & antidiffusion_bounds,
  const Vector<double> & cumulative_antidiffusion_vector,
  SparseMatrix<double> & limiter_matrix)
{
  // pass index
  unsigned int pass_index = 1;

  // infinite loop
  while (true)
  {
    // determine if and how much antidiffusion was accepted in this pass
    const double antidiffusion_added = compute_added_antidiffusion();

    bool no_more_antidiffusion;

    // break if no antidiffusion was accepted in this pass
    if (no_more_antidiffusion)
      break;

    // increment pass index
    pass_index++;
  }
}
