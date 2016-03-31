/**
 * \file Limiter.cc
 * \brief Provides the function definitions for the Limiter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 * \param[in] report_antidiffusion_  flag to report amount of accepted
 *            antidiffusion
 */
template <int dim>
Limiter<dim>::Limiter(const unsigned int & n_dofs_,
                      const bool & report_antidiffusion_)
  : n_dofs(n_dofs_),
    zero_vector(n_dofs_),
    report_antidiffusion(report_antidiffusion_)
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

  // report antidiffusion if requested
  if (report_antidiffusion)
  {
    // compute percentage of possible antidiffusion that was accepted
    const double total_possible_antidiffusion =
      compute_total_possible_antidiffusion(antidiffusion_matrix);
    const double total_antidiffusion =
      compute_total_antidiffusion(antidiffusion_matrix, limiter_matrix);
    const double percentage_antidiffusion =
      total_antidiffusion / total_possible_antidiffusion * 100.0;

    printf("      %.2f%% of antidiffusion accepted.\n", percentage_antidiffusion);
  }
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

/**
 * \brief Computes the total possible antidiffusion, used for reporting.
 *
 * Considering that antidiffusive fluxes come in pairs of equal and opposite
 * values, the total possible antidiffusion is defined to be
 * \f[
 *   p^{total} = \frac{1}{2}\sum\limits_{i,j}|P_{i,j}| \,.
 * \f]
 *
 * \param[in] antidiffusion_matrix  matrix of antidiffusion fluxes
 *            \f$\mathbf{P}\f$
 */
template <int dim>
double Limiter<dim>::compute_total_possible_antidiffusion(
  const SparseMatrix<double> & antidiffusion_matrix) const
{
  // initialize sum to zero
  double total_antidiffusion = 0.0;

  // loop over rows
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // get first and end iterator for row
    SparseMatrix<double>::const_iterator it = antidiffusion_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_end = antidiffusion_matrix.end(i);

    // loop over row entries
    for (; it != it_end; ++it)
      total_antidiffusion += std::abs(it->value());
  }

  // divide by 2 since each flux magnitude was added twice
  total_antidiffusion *= 0.5;

  return total_antidiffusion;
}

/**
 * \brief Computes the total antidiffusion accepted.
 *
 * Considering that antidiffusive fluxes come in pairs of equal and opposite
 * values, the total antidiffusion is defined to be
 * \f[
 *   \bar{p}^{total} = \frac{1}{2}\sum\limits_{i,j}|L_{i,j}P_{i,j}| \,.
 * \f]
 *
 * \param[in] antidiffusion_matrix  matrix of antidiffusion fluxes
 *            \f$\mathbf{P}\f$
 * \param[in] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$
 */
template <int dim>
double Limiter<dim>::compute_total_antidiffusion(
  const SparseMatrix<double> & antidiffusion_matrix,
  const SparseMatrix<double> & limiter_matrix) const
{
  // initialize sum to zero
  double total_antidiffusion = 0.0;

  // loop over rows
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // get first and end iterator for row
    SparseMatrix<double>::const_iterator it_flux = antidiffusion_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_end = antidiffusion_matrix.end(i);
    SparseMatrix<double>::const_iterator it_lim = limiter_matrix.begin(i);

    // loop over row entries
    for (; it_flux != it_end; ++it_flux, ++it_lim)
      total_antidiffusion += std::abs(it_lim->value() * it_flux->value());
  }

  // divide by 2 since each flux magnitude was added twice
  total_antidiffusion *= 0.5;

  return total_antidiffusion;
}
