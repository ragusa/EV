/**
 * \file Limiter.cc
 * \brief Provides the function definitions for the Limiter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 */
template <int dim>
Limiter<dim>::Limiter(const unsigned int & n_dofs_)
  : n_dofs(n_dofs_), zero_vector(n_dofs_)
{
}

/**
 * \brief Computes the limiting coefficient matrix.
 *
 * This function calls the more general function defined in derived classes,
 * which allows a cumulative antidiffusion vector to be passed. This function
 * passes a zero vector for this parameter.
 *
 * \param[in] antidiffusion_matrix  matrix of antidiffusion fluxes
 *            \f$\mathbf{P}\f$
 * \param[in] antidiffusion_bounds  lower and upper bounds of antidiffusive
 *            fluxes into each node \f$\mathbf{Q}^-\f$ and \f$\mathbf{Q}^+\f$
 * \param[out] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$
 */
template <int dim>
void Limiter<dim>::compute_limiter_matrix(
  const SparseMatrix<double> & antidiffusion_matrix,
  const DoFBounds<dim> & antidiffusion_bounds,
  SparseMatrix<double> & limiter_matrix)
{
  // call function in derived class
  compute_limiter_matrix(
    antidiffusion_matrix, antidiffusion_bounds, zero_vector, limiter_matrix);
}

/**
 * \brief Applies limiting coefficients to an antidiffusion matrix.
 *
 * Limiting coefficients are applied via an element-wise multiplication of the
 * limiter matrix \f$\mathbf{L}\f$ and the antidiffusion matrix \f$\mathbf{P}\f$:
 * \f[
 *   P_{i,j} \gets L_{i,j}P_{i,j} \,.
 * \f]
 *
 * \pre This function assumes that \f$\mathbf{L}\f$ and \f$\mathbf{P}\f$ have
 *      the same sparsity pattern.
 *
 * \param[in] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$
 * \param[inout] antidiffusion_matrix  matrix of antidiffusion fluxes
 *               \f$\mathbf{P}\f$
 */
template <int dim>
void Limiter<dim>::apply_limiter_matrix(
  const SparseMatrix<double> & limiter_matrix,
  SparseMatrix<double> & antidiffusion_matrix) const
{
  // loop over rows
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // get first and end iterator for row
    SparseMatrix<double>::const_iterator it_lim = limiter_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_end = limiter_matrix.end(i);
    SparseMatrix<double>::iterator it_flux = antidiffusion_matrix.begin(i);

    // loop over row entries
    for (; it_lim != it_end; ++it_lim, ++it_flux)
      it_flux->value() = it_lim->value() * it_flux->value();
  }
}
