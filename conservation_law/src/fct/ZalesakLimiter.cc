/**
 * \file ZalesakLimiter.cc
 * \brief Provides the function definitions for the ZalesakLimiter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 */
ZalesakLimiter::ZalesakLimiter(const unsigned int & n_dofs_) : Limiter(n_dofs_)
{
  // resize limiter vectors
  negative_limiter_vector.reinit(n_dofs)
  positive_limiter_vector.reinit(n_dofs)
}

/**
 * \brief Computes the limiting coefficient matrix.
 *
 * The Zalesak limiter limits antidiffusive fluxes \f$P_{i,j}\f$ to satisfy the
 * following conditions:
 * \f[
 *   Q_i^- \le L_{i,j}P_{i,j} \le Q_i^+ \quad \forall i \,,
 * \f]
 * where \f$Q_i^-\f$ and \f$Q_i^+\f$ are antidiffusion bounds for degree of
 * freedom \f$i\f$, and \f$L_{i,j}\f$ is the limiting coefficient corresponding
 * to the antidiffusive flux \f$P_{i,j}\f$.
 *
 * \note This function assumes that if a degree of freedom has Dirichlet
 *       constraints, then the antidiffusion input vectors \f$\mathbf{Q}^-\f$
 *       and \f$\mathbf{Q}^+\f$ have been treated accordingly so as to
 *       prevent limitation of antidiffusion fluxes to and from these
 *       nodes; this amounts to setting \f$Q_i^-=-c_{large}\f$ and
 *       \f$Q_i^+=c_{large}\f$ for all Dirichlet nodes \f$i\f$, where
 *       \f$c_{large}\f$ is a positive number with a large magnitude.
 *
 * \param[in] antidiffusion_matrix  matrix of antidiffusion fluxes
 *            \f$\mathbf{P}\f$ 
 * \param[in] antidiffusion_bounds  lower and upper bounds of antidiffusive
 *            fluxes into each node \f$\mathbf{Q}^-\f$ and \f$\mathbf{Q}^+\f$
 * \param[out] limiter_matrix  matrix of limiting coeffients \f$\mathbf{L}\f$ 
 */
void ZalesakLimiter::compute_limiter(const SparseMatrix<double> & antidiffusion_matrix,
  const DoFBounds & antidiffusion_bounds,
  SparseMatrix<double> & limiter_matrix) const
{
  // reset limiter matrix
  limiter_matrix = 0;

  // compute antidiffusion limiter vectors L- and L+
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    // get first and end iterator for row
    SparseMatrix<double>::const_iterator it = antidiffusion_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_end = antidiffusion_matrix.end(i);

    // loop over row entries to compute negative and positive antidiffusive
    // flux sums P- and P+
    double P_negative = 0.0;
    double P_positive = 0.0;
    for (; it != it_end; ++it)
    {
      const double Pij = it->value();
      P_negative += std::min(0.0, Pij);
      P_positive += std::max(0.0, Pij);
    }

    // compute L-[i]
    if (P_negative != 0.0)
      negative_limiter_vector[i] = std::min(1.0, antidiffusion_bounds.get_min(i) / P_negative);
    else
      negative_limiter_vector[i] = 1.0;

    // compute L+[i]
    if (P_positive != 0.0)
      positive_limiter_vector[i] = std::min(1.0, antidiffusion_bounds.get_max(i) / P_positive);
    else
      positive_limiter_vector[i] = 1.0;
  }

  // compute limited flux correction sum
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    // get first and end iterator for row
    SparseMatrix<double>::const_iterator it_flux = antidiffusion_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_end = antidiffusion_matrix.end(i);
    SparseMatrix<double>::iterator it_lim = limiter_matrix.begin(i);

    // loop over row entries
    for (; it_flux != it_end; ++it_flux, ++it_lim)
    {
      // get column index
      const unsigned int j = it_flux->column();

      // get antidiffusion flux value
      const double Pij = it_flux->value();

      // compute limiting coefficient
      double Lij;
      if (Pij >= 0.0)
        Lij = std::min(positive_limiter_vector[i], negative_limiter_vector[j]);
      else
        Lij = std::min(positive_limiter_vector[j], negative_limiter_vector[i]);

      // store limiting coefficient
      it_lim->value() = Lij;
    }
  }
}
