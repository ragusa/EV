/**
 * \file MultipassLimiter.cc
 * \brief Provides the function definitions for the MultipassLimiter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 * \param[in] limiter_  limiter for which multiple passes are to be made
 * \param[in] sparsity_pattern_  sparsity pattern
 * \param[in] percent_tolerance_  percent tolerance for ending passes
 * \param[in] report_antidiffusion_  flag to report amount of accepted
 *            antidiffusion
 */
template <int dim>
MultipassLimiter<dim>::MultipassLimiter(
  const unsigned int & n_dofs_,
  std::shared_ptr<Limiter<dim>> limiter_,
  const std::shared_ptr<SparsityPattern> sparsity_pattern_,
  const double & percent_tolerance_,
  const bool & report_antidiffusion_)
  : Limiter<dim>(n_dofs_, report_antidiffusion_),
    limiter(limiter_),
    sparsity_pattern(sparsity_pattern_),
    percent_tolerance(percent_tolerance_),
    accepted_antidiffusion(n_dofs_)
{
  // initialize matrices with sparsity pattern
  remainder_matrix.reinit(*sparsity_pattern);
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
 *            iterations \f$\bar{\mathbf{p}}^{(\ell-1)}\f$
 * \param[out] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$
 */
template <int dim>
void MultipassLimiter<dim>::compute_limiter_matrix(
  const SparseMatrix<double> & antidiffusion_matrix,
  const DoFBounds<dim> & antidiffusion_bounds,
  const Vector<double> &,
  SparseMatrix<double> & limiter_matrix)
{
  // pass index
  unsigned int pass_index = 1;

  // initialize remainder matrix and accepted antidiffusion vector
  remainder_matrix.copy_from(antidiffusion_matrix);
  accepted_antidiffusion = 0;

  // compute total possible antidiffusion; used for reporting
  const double total_possible_antidiffusion =
    this->compute_total_possible_antidiffusion(antidiffusion_matrix);

  // infinite loop
  while (true)
  {
    // call limiter
    limiter->compute_limiter_matrix(remainder_matrix,
                                    antidiffusion_bounds,
                                    accepted_antidiffusion,
                                    limiter_matrix);

    // determine if and how much antidiffusion was accepted in this pass
    const double antidiffusion_added =
      this->compute_total_antidiffusion(remainder_matrix, limiter_matrix);
    const double percentage_antidiffusion =
      antidiffusion_added / total_possible_antidiffusion * 100.0;

    // report antidiffusion added in this iteration
    if (this->report_antidiffusion)
      printf("      Pass %i: %.2f%% of antidiffusion accepted.\n",
             pass_index,
             percentage_antidiffusion);

    // break if no antidiffusion was accepted in this pass
    if (percentage_antidiffusion < percent_tolerance)
      break;

    // update remainder antidiffusion matrix
    update_remainder_matrix_and_accepted_antidiffusion(limiter_matrix);

    // increment pass index
    pass_index++;
  }

  // compute total limiting coefficients for all passes
  compute_total_limiter_matrix(
    antidiffusion_matrix, remainder_matrix, limiter_matrix);
}

/**
 * \brief Updates the remainder antidiffusion matrix and the accepted
 *        antidiffusion vector.
 *
 * The update to each remainder antidiffusion matrix entry is performed as
 * follows:
 * \f[
 *   \Delta P_{i,j}^{(\ell+1)} = (1 - L_{i,j}^{(\ell)})\Delta P_{i,j}^{(\ell)} \,.
 * \f]
 * The update to each accepted antidiffusion vector entry is performed as
 * follows:
 * \f[
 *   \bar{p}_i^{(\ell+1)} = \bar{p}_i^{(\ell)}
 *     + \sum\limits_j L_{i,j}^{(\ell)}\Delta P_{i,j}^{(\ell)} \,.
 * \f]
 *
 * \param[in] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$
 */
template <int dim>
void MultipassLimiter<dim>::update_remainder_matrix_and_accepted_antidiffusion(
  const SparseMatrix<double> & limiter_matrix)
{
  // loop over degrees of freedom
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    // create sparse matrix iterators
    SparseMatrix<double>::iterator it_remainder = remainder_matrix.begin(i);
    SparseMatrix<double>::iterator it_end = remainder_matrix.end(i);
    SparseMatrix<double>::const_iterator it_lim = limiter_matrix.begin(i);

    // loop over entries in row
    for (; it_remainder != it_end; ++it_remainder, ++it_lim)
    {
      const double DPij = it_remainder->value();
      const double Lij = it_lim->value();

      // update remainder antidiffusion matrix entry
      it_remainder->value() *= (1.0 - Lij);

      // update accepted antidiffusion
      accepted_antidiffusion[i] += Lij * DPij;
    }
  }
}

/**
 * \brief Computes the total limiting coefficients for all passes.
 *
 * This function uses the remainder antidiffusive fluxes after the last pass to
 * compute total limiting coefficients:
 * \f[
 *   \Delta P_{i,j} = P_{i,j} - L_{i,j}^{total}P_{i,j} \,,
 * \f]
 * and thus the total limiting coefficients are
 * \f[
 *   L_{i,j}^{total} = 1 - \frac{\Delta P_{i,j}}{P_{i,j}} \,.
 * \f]
 * In case there would be division by zero (when \f$P_{i,j} = 0\f$),
 * \f$L_{i,j}^{total}\f$ is taken to be one.
 *
 * \param[in] antidiffusion_matrix  matrix of antidiffusion fluxes
 *            \f$\mathbf{P}\f$
 * \param[in] remainder_antidiffusion_matrix  matrix of remainder antidiffusion
 *            fluxes \f$\Delta\mathbf{P}\f$
 * \param[out] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$
 */
template <int dim>
void MultipassLimiter<dim>::compute_total_limiter_matrix(
  const SparseMatrix<double> & antidiffusion_matrix,
  const SparseMatrix<double> & remainder_antidiffusion_matrix,
  SparseMatrix<double> & limiter_matrix) const
{
  // loop over degrees of freedom
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    // create sparse matrix iterators
    SparseMatrix<double>::iterator it_lim = limiter_matrix.begin(i);
    SparseMatrix<double>::iterator it_end = limiter_matrix.end(i);
    SparseMatrix<double>::const_iterator it_flux = antidiffusion_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_remainder =
      remainder_antidiffusion_matrix.begin(i);

    // loop over entries in row
    for (; it_lim != it_end; ++it_lim, ++it_flux, ++it_remainder)
    {
      const double Pij = it_flux->value();
      const double DPij = it_remainder->value();
      if (Pij == 0.0)
        it_lim->value() = 1.0;
      else
        it_lim->value() = 1.0 - DPij / Pij;
    }
  }
}
