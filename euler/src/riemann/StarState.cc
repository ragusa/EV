/**
 * \file StarState.cc
 * \brief Provides the function definitions for the StarState
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] max_wave_speed_ max wave speed object
 * \param[in] triangulation_ triangulation
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] n_components_ number of solution components
 * \param[in] n_dofs number of degrees of freedom in vector system
 */
template <int dim>
StarState<dim>::StarState(
  const std::shared_ptr<MaxWaveSpeed<dim>> & max_wave_speed_,
  const Triangulation<dim> & triangulation_,
  const QGauss<dim> & cell_quadrature_,
  const unsigned int & n_components_,
  const unsigned int & n_dofs_)
  : ArtificialDiffusion<dim>(),
    max_wave_speed(max_wave_speed_),
    fe_scalar(1),
    dof_handler_scalar(triangulation_),
    cell_quadrature(&cell_quadrature_),
    n_components(n_components_),
    n_dofs(n_dofs_),
    n_dofs_scalar(n_dofs / n_components),
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
void StarState<dim>::reinitialize()
{
  // distribute degrees of freedom for scalar case
  dof_handler_scalar.clear();
  dof_handler_scalar.distribute_dofs(fe_scalar);

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
 * \brief Computes artificial diffusion matrix.
 *
 * This matrix is computed as
 * \f[
 *   D_{i,j} = \left\{\begin{array}{c l}
 *     -\max\left(
 *     \lambda_{max}(\mathbf{n}_{i,j},\mathbf{u}_i^n,\mathbf{u}_j^n)
 *     \|\mathbf{c}_{i,j}\|_{\ell^2},
 *     \lambda_{max}(\mathbf{n}_{j,i},\mathbf{u}_j^n,\mathbf{u}_i^n)
 *     \|\mathbf{c}_{j,i}\|_{\ell^2}
 *     \right) & j\ne i \\
 *     -\sum\limits_{k\ne i}D_{i,k} & j = i
 *     \end{array}\right. \,.
 * \f]
 *
 * \param[in] solution solution vector
 * \param[in] viscosity pointer to viscosity cell map
 * \param[out] diffusion_matrix diffusion matrix
 */
template <int dim>
void StarState<dim>::compute_diffusion_matrix(
  const Vector<double> & solution,
  const std::shared_ptr<Viscosity<dim>>,
  SparseMatrix<double> & diffusion_matrix)
{
  // reset diffusion matrix to zero
  diffusion_matrix = 0;

  // loop over rows of a scalar matrix
  for (unsigned int i = 0; i < n_dofs_scalar; ++i)
  {
    // NOTE: for the vector matrix, it might be advantageous to use iterators
    // to access the entries; iterators can be incremented manually as each
    // component is set

    // intialize iterators
    std::vector<SparseMatrix<double>::const_iterator> it;
    for (unsigned int d = 0; d < dim; ++d)
      it.push_back(gradients[d].begin(i));
    SparseMatrix<double>::const_iterator it_end = gradients[0].end(i);
    SparseMatrix<double>::const_iterator it_norm = gradient_norms.begin(i);

    // loop over row entries of scalar matrices
    for (; it[0] != it_end; ++it_norm)
    {
      // get row and column indices
      const unsigned int i = it[0]->row();
      const unsigned int j = it[0]->column();

      // compute indices of first component in vector system
      const unsigned int i0 = i * n_components;
      const unsigned int j0 = j * n_components;

      // make sure this is an off-diagonal entry; diagonal computed afterward
      if (j != i)
      {
        // NOTE: Is it safe to assume that after resetting the matrix the zero
        // at the beginning, all entries are exactly zero and thus the floating
        // point comparison that follows is valid? If not, perhaps a boolean
        // matrix should be employed to make this check.

        // check if entry has already been computed
        if (diffusion_matrix(i0, j0) == 0.0)
        {
          // get solution at each support point
          std::vector<double> solution_i(n_components);
          std::vector<double> solution_j(n_components);
          for (unsigned int m = 0; m < n_components; ++m)
          {
            solution_i[m] = solution[i * n_components + m];
            solution_j[m] = solution[j * n_components + m];
          }

          // NOTE: Here I should check if either i or j corresponds to an interior
          // support point. If so, then as given in Remark 3.2 of Guermond's
          // invariant domain paper, lambda_max(i,j) = lambda_max(j,i), so there
          // is
          // no need to compute both max wave speeds here.

          // create tensor for normal vector
          Tensor<1, dim> normal_ij;
          Tensor<1, dim> normal_ji;
          for (unsigned int d = 0; d < dim; ++d)
          {
            // note that gradients matrix has been normalized already
            normal_ij[d] = it[d]->value();
            normal_ji[d] = gradients[d](j, i);
          }

          // compute maximum wave speed
          const double max_speed_ij =
            max_wave_speed->compute(solution_i, solution_j, normal_ij);
          const double max_speed_ji =
            max_wave_speed->compute(solution_j, solution_i, normal_ji);

          // get gradient norms
          const double c_ij = it_norm->value();
          const double c_ji = gradient_norms(j, i);

          // compute diffusion entry value (note D(i,j) = D(j,i))
          const double D_ij = -std::max(max_speed_ij * c_ij, max_speed_ji * c_ji);

          // loop over solution components (all share same diffusion matrix)
          for (unsigned int m = 0; m < n_components; ++m)
          {
            // compute indices for component
            const unsigned int im = i * n_components + m;
            const unsigned int jm = j * n_components + m;

            // set diffusion matrix entries
            diffusion_matrix.set(im, jm, D_ij);
            diffusion_matrix.set(jm, im, D_ij);

            // subtract contribution from diagonals
            diffusion_matrix.add(im, im, -D_ij);
            diffusion_matrix.add(jm, jm, -D_ij);
          }
        }
      } // if j != i

      // increment iterators
      for (unsigned int d = 0; d < dim; ++d)
        ++it[d];
    } // end row iterator loop
  }   // end loop over rows
}

/**
 * \brief Computes gradients and normal vectors along each edge.
 */
template <int dim>
void StarState<dim>::compute_gradients_and_normals()
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
