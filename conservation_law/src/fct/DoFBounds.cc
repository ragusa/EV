/**
 * \file DoFBounds.cc
 * \brief Provides the function definitions for the DoFBounds class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_dofs_  number of degrees of freedom
 */
DoFBounds<dim>::DoFBounds(const unsigned int & n_dofs_) : n_dofs(n_dofs_)
{
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
    lower_bound[i] = -1.0e15;
    upper_bound[i] = 1.0e15;
  }
}
