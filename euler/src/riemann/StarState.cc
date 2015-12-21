/**
 * \file StarState.cc
 * \brief Provides the function definitions for the StarState
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] gradient_matrix_ pointer to gradient matrix object
 * \param[in] dof_handler_ degree of freedom handler for vector case
 * \param[in] triangulation_ triangulation
 * \param[in] n_components_ number of solution components
 */
template <int dim>
StarState<dim>::StarState(
  const std::shared_ptr<GradientMatrix<dim>> & gradient_matrix_,
  const DoFHandler<dim> & dof_handler_,
  const Triangulation<dim> & triangulation_,
  const unsigned int & n_components_)
  : gradient_matrix(gradient_matrix_),
    fe_scalar(1),
    dof_handler(&dof_handler_),
    dof_handler_scalar(triangulation_),
    n_components(n_components_)
{
}

/**
 * \brief Reinitializes.
 */
template <int dim>
void StarState<dim>::reinitialize()
{
  // distribute degrees of freedom for scalar case
  dof_handler_scalar.clear();
  dof_handler_scalar.distribute_dofs(fe_scalar);
  n_dofs_scalar = dof_handler_scalar.n_dofs();

  // create sparsity pattern for a scalar matrix
  DynamicSparsityPattern dsp_scalar(n_dofs_scalar);
  DoFTools::make_sparsity_pattern(dof_handler_scalar, dsp_scalar);
  sparsity_scalar.copy_from(dsp_scalar);

  // number of degrees of freedom for vector case
  n_dofs = dof_handler->n_dofs();

  // create sparity pattern for a vector matrix and initialize matrix
  // for star states
  DynamicSparsityPattern dsp(n_dofs);
  DoFTools::make_sparsity_pattern(*dof_handler, dsp);
  sparsity.copy_from(dsp);
  star_state_matrix.reinit(sparsity);
}

/**
 * \brief Computes star states for each pair of nodal solutions
 *        \f$\mathbf{u}_i\f$ and \f$\mathbf{u}_j\f$ with shared support.
 *
 * \param[in] solution solution vector
 */
template <int dim>
void StarState<dim>::compute_star_states(const Vector<double> & solution)
{
  // loop over entries of scalar matrix
  for (unsigned int i = 0; i < n_dofs_scalar; ++i)
  {
    // get solution at each support point
    std::vector<double> solution_i(n_components);
    for (unsigned int m = 0; m < n_components; ++m)
      solution_i[m] = solution[i * n_components + m];

    // iterate over scalar sparsity pattern
    SparsityPattern::iterator it = sparsity_scalar.begin(i);
    SparsityPattern::iterator it_end = sparsity_scalar.end(i);
    for (; it != it_end; ++it)
    {
      // get column
      const unsigned int j = it->column();

      // get solution at each support point
      std::vector<double> solution_j(n_components);
      for (unsigned int m = 0; m < n_components; ++m)
        solution_j[m] = solution[j * n_components + m];

      // get normal for i and j
      auto normal_ij = gradient_matrix->get_normal(i, j);

      // compute star state solution
      auto solution_star = compute(solution_i, solution_j, normal_ij);

      // store star state solution in vector matrix
      for (unsigned int m = 0; m < n_components; ++m)
      {
        // compute indices
        const unsigned int i_m = i * n_components + m;
        const unsigned int j_m = j * n_components + m;

        // store value
        star_state_matrix.set(i_m, j_m, solution_star[m]);
      }
    }
  }
}
