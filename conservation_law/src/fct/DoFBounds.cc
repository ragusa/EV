/**
 * \file DoFBounds.cc
 * \brief Provides the function definitions for the DoFBounds class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 */
DoFBounds::DoFBounds(const unsigned int & n_dofs_) : n_dofs(n_dofs_)
{
  // resize bounds vectors
  lower.reinit(n_dofs);
  upper.reinit(n_dofs);
}

/**
 * \brief Widens bounds \f$\mathbf{W}^\pm\f$ to another set of bounds
 *        \f$\mathbf{X}^\pm\f$.
 *
 * This function replaces the upper bound by the maximum of the upper bound
 * and the other upper bound:
 * \f[
 *   W_i^+ \gets \max(W_i^+, X_i^+) \,,
 * \f]
 * and similarly, replaces the lower bound by the minimum of the lower bound
 * and the other lower bound:
 * \f[
 *   W_i^- \gets \min(W_i^-, X_i^-) \,.
 * \f]
 *
 * \param[in] other_dof_bounds  other degree of freedom bounds
 *            \f$\mathbf{X}^\pm\f$
 */
void DoFBounds::widen(const DoFBounds & other_dof_bounds)
{
  // loop over degrees of freedom
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    lower[i] = std::min(lower[i], other_dof_bounds.lower[i]);
    upper[i] = std::max(upper[i], other_dof_bounds.upper[i]);
  }
}

/**
 * \brief Shrinks bounds \f$\mathbf{W}^\pm\f$ to another set of bounds
 *        \f$\mathbf{X}^\pm\f$.
 *
 * This function replaces the upper bound by the minimum of the upper bound
 * and the other upper bound:
 * \f[
 *   W_i^+ \gets \min(W_i^+, Y_i^+) \,,
 * \f]
 * and similarly, replaces the lower bound by the maximum of the lower bound
 * and the other lower bound:
 * \f[
 *   W_i^- \gets \max(W_i^-, Y_i^-) \,.
 * \f]
 *
 * \param[in] other_dof_bounds  other degree of freedom bounds
 *            \f$\mathbf{X}^\pm\f$
 */
void DoFBounds::shrink(const DoFBounds & other_dof_bounds)
{
  // loop over degrees of freedom
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    lower[i] = std::max(lower[i], other_dof_bounds.lower[i]);
    upper[i] = std::min(upper[i], other_dof_bounds.upper[i]);
  }
}

/**
 * \brief Removes bounds for a list of degrees of freedom.
 *
 * This function removes bounds for a degree of freedom by making the lower
 * bound a very large negative number and the upper bound a very large
 * positive number.
 *
 * \param[in] dof_indices  vector of degree of freedom indices for which bounds
 *            are to be removed
 */
void DoFBounds::remove_bounds(const std::vector<double> & dof_indices)
{
  // get number of indices for which bounds are to be removed
  const unsigned int n_indices = dof_indices.size();

  // loop over indices for which bounds are to be removed
  for (unsigned int j = 0; j < n_indices; ++j)
  {
    // get degree of freedom index
    const unsigned int i = dof_indices[j];

    // make bounds very large to effectively remove them
    lower[i] = -1.0e15;
    upper[i] = 1.0e15;
  }
}
