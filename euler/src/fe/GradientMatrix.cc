/**
 * \file GradientMatrix.cc
 * \brief Provides the function definitions for the GradientMatrix
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] triangulation_ triangulation
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] n_components_ number of solution components
 */
template <int dim>
GradientMatrix<dim>::GradientMatrix(const Triangulation<dim> & triangulation_,
                                    const QGauss<dim> & cell_quadrature_,
                                    const unsigned int & n_components_)
  : fe_scalar(1),
    dof_handler_scalar(triangulation_),
    cell_quadrature(&cell_quadrature_),
    n_components(n_components_),
    dofs_per_cell_scalar(fe_scalar.dofs_per_cell),
    n_q_points_cell(cell_quadrature_.size())
{
  // reinitialize
  reinitialize();
}

/**
 * \brief Reinitializes.
 */
template <int dim>
void GradientMatrix<dim>::reinitialize()
{
  // distribute degrees of freedom for scalar case
  dof_handler_scalar.clear();
  dof_handler_scalar.distribute_dofs(fe_scalar);

  // get number of degrees of freedom
  n_dofs_scalar = dof_handler_scalar.n_dofs();
  n_dofs = n_dofs_scalar * n_components;

  // create sparsity pattern and reinitialize sparse matrices
  DynamicSparsityPattern dsp(n_dofs_scalar);
  DoFTools::make_sparsity_pattern(dof_handler_scalar, dsp);
  sparsity.copy_from(dsp);
  for (unsigned int d = 0; d < dim; ++d)
    gradients[d].reinit(sparsity);
  gradient_norms.reinit(sparsity);

  // compute gradients
  compute_gradients_and_normals();
}

/**
 * \brief Computes gradients and normal vectors along each edge.
 */
template <int dim>
void GradientMatrix<dim>::compute_gradients_and_normals()
{
  // reset matrices
  for (unsigned int d = 0; d < dim; ++d)
    gradients[d] = 0;
  gradient_norms = 0;

  FEValues<dim> fe_values(fe_scalar,
                          *cell_quadrature,
                          update_values | update_gradients | update_JxW_values);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell_scalar);

  // loop over cells
  Cell cell = dof_handler_scalar.begin_active(), endc = dof_handler_scalar.end();
  for (; cell != endc; ++cell)
  {
    // reinitialize fe values for cell
    fe_values.reinit(cell);

    // get DoF indices
    cell->get_dof_indices(local_dof_indices);

    // loop over dimension
    for (unsigned int d = 0; d < dim; ++d)
      // loop over test function i
      for (unsigned int i = 0; i < dofs_per_cell_scalar; ++i)
        // loop over test function j
        for (unsigned int j = 0; j < dofs_per_cell_scalar; ++j)
        {
          // initialize integral to zero
          double entry_cell_i_j = 0.0;

          // loop over quadrature points
          for (unsigned int q = 0; q < n_q_points_cell; ++q)
            entry_cell_i_j += fe_values.shape_value(i, q) *
              fe_values.shape_grad(j, q)[d] * fe_values.JxW(q);

          // add contribution to global matrix
          gradients[d].add(
            local_dof_indices[i], local_dof_indices[j], entry_cell_i_j);
        }
  }

  // initialize matrix iterators
  SparseMatrix<double>::iterator it = gradient_norms.begin();
  SparseMatrix<double>::iterator it_end = gradient_norms.end();
  std::vector<SparseMatrix<double>::iterator> it_grad;
  for (unsigned int d = 0; d < dim; ++d)
    it_grad.push_back(gradients[d].begin());

  for (; it != it_end; ++it)
  {
    // compute gradient norms
    for (unsigned int d = 0; d < dim; ++d)
      it->value() += it_grad[d]->value() * it_grad[d]->value();
    it->value() = std::sqrt(it->value());

    // compute normalized gradient vector
    for (unsigned int d = 0; d < dim; ++d)
    {
      // normalize this component of gradient
      it_grad[d]->value() = it_grad[d]->value() / it->value();

      // increment iterator
      ++it_grad[d];
    }
  }
}

/**
 * \brief Gets the normal vector for degrees of freedom i and j.
 *
 * \param[in] i row index
 * \param[in] j column index
 *
 * \return normal vector for degrees of freedom i and j
 */
template <int dim>
Tensor<1, dim> GradientMatrix<dim>::get_normal(const unsigned int & i,
                                               const unsigned int & j) const
{
  Tensor<1, dim> normal;
  for (unsigned int d = 0; d < dim; ++d)
    normal[d] = gradients[d](i, j);

  return normal;
}

/**
 * \brief Gets the normal vectors for row i
 *
 * \param[in] i row index
 *
 * \return normal vectors for row i
 */
template <int dim>
std::vector<Tensor<1, dim>> GradientMatrix<dim>::get_normals(
  const unsigned int & i) const
{
  std::vector<Tensor<1, dim>> normals;

  // initialize iterators
  std::vector<SparseMatrix<double>::const_iterator> it;
  for (unsigned int d = 0; d < dim; ++d)
    it.push_back(gradients[d].begin(i));
  SparseMatrix<double>::const_iterator it_end = gradients[0].end(i);

  // loop over columns of row
  for (; it[0] != it_end;)
  {
    // get normal
    Tensor<1, dim> normal;
    for (unsigned int d = 0; d < dim; ++d)
      normal[d] = it[d]->value();

    // add normal to vector
    normals.push_back(normal);

    // increment iterators
    for (unsigned int d = 0; d < dim; ++d)
      ++it[d];
  }

  return normals;
}

/**
 * \brief Gets the norm of the gradient vector for degrees of freedom i and j.
 *
 * \param[in] i row index
 * \param[in] j column index
 *
 * \return norm of the gradient vector for degrees of freedom i and j
 */
template <int dim>
double GradientMatrix<dim>::get_gradient_norm(const unsigned int & i,
                                              const unsigned int & j) const
{
  return gradient_norms(i, j);
}

/**
 * \brief Gets the norms of the gradient vectors for row i.
 *
 * \param[in] i row index
 *
 * \return norms of the gradient vectors for row i
 */
template <int dim>
std::vector<double> GradientMatrix<dim>::get_gradient_norms(
  const unsigned int & i) const
{
  std::vector<double> norms;

  SparseMatrix<double>::const_iterator it = gradient_norms.begin(i);
  SparseMatrix<double>::const_iterator it_end = gradient_norms.end(i);
  for (; it != it_end; ++it)
    norms.push_back(it->value());

  return norms;
}

/**
 * \brief Gets the nonzero column indices for row i.
 *
 * \param[in] i row index
 *
 * \return nonzero column indices for row i
 */
template <int dim>
std::vector<unsigned int> GradientMatrix<dim>::get_column_indices(
  const unsigned int & i) const
{
  std::vector<unsigned int> column_indices;

  SparseMatrix<double>::const_iterator it = gradient_norms.begin(i);
  SparseMatrix<double>::const_iterator it_end = gradient_norms.end(i);
  for (; it != it_end; ++it)
    column_indices.push_back(it->column());

  return column_indices;
}
