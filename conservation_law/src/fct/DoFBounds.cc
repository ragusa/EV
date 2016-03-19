/**
 * \file DoFBounds.cc
 * \brief Provides the function definitions for the DoFBounds class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 */
template <int dim>
DoFBounds<dim>::DoFBounds(const DoFHandler<dim> & dof_handler_,
                          const FESystem<dim> & fe_)
  : dof_handler(&dof_handler_),
    fe(&fe_),
    n_dofs(dof_handler_.n_dofs()),
    dofs_per_cell(fe_.dofs_per_cell)
{
  // resize bounds vectors
  lower.reinit(n_dofs);
  upper.reinit(n_dofs);
}

/**
 * \brief Computes the minimum and maximum of a DoF vector
 *        in the neighborhood of each degree of freedom.
 *
 * \param[in] dof_vector       DoF vector at which to evaluate min/max
 * \param[out] min_dof_vector  min of DoF vector in neighborhood of each DoF
 * \param[out] max_dof_vector  max of DoF vector in neighborhood of each DoF
 */
template <int dim>
void DoFBounds<dim>::compute_min_max_dof_vector(
  const Vector<double> & dof_vector,
  Vector<double> & min_dof_vector,
  Vector<double> & max_dof_vector) const
{
  // initialize min and max
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    min_dof_vector(i) = dof_vector(i);
    max_dof_vector(i) = dof_vector(i);
  }

  // loop over cells
  Cell cell = dof_handler->begin_active(), endc = dof_handler->end();
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // find min and max values on cell
    double max_cell = dof_vector(local_dof_indices[0]);
    double min_cell = dof_vector(local_dof_indices[0]);
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      double value_j = dof_vector(local_dof_indices[j]);
      max_cell = std::max(max_cell, value_j);
      min_cell = std::min(min_cell, value_j);
    }

    // update the max and min values of neighborhood of each dof
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      unsigned int i = local_dof_indices[j]; // global index
      min_dof_vector(i) = std::min(min_dof_vector(i), min_cell);
      max_dof_vector(i) = std::max(max_dof_vector(i), max_cell);
    }
  }
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
template <int dim>
void DoFBounds<dim>::widen(const DoFBounds & other_dof_bounds)
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
template <int dim>
void DoFBounds<dim>::shrink(const DoFBounds & other_dof_bounds)
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
template <int dim>
void DoFBounds<dim>::remove_bounds(const std::vector<double> & dof_indices)
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
