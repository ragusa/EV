/**
 * \file DomainInvariantDiffusion.cc
 * \brief Provides the function definitions for the DomainInvariantDiffusion
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] max_wave_speed_ max wave speed object
 * \param[in] gradient_matrix_ gradient matrix
 * \param[in] n_components_ number of solution components
 * \param[in] n_dofs number of degrees of freedom in vector system
 */
template <int dim>
DomainInvariantDiffusion<dim>::DomainInvariantDiffusion(
  const std::shared_ptr<MaxWaveSpeed<dim>> & max_wave_speed_,
  const std::shared_ptr<GradientMatrix<dim>> & gradient_matrix_,
  const unsigned int & n_components_,
  const unsigned int & n_dofs_)
  : ArtificialDiffusion<dim>(),
    max_wave_speed(max_wave_speed_),
    gradient_matrix(gradient_matrix_),
    n_components(n_components_),
    n_dofs_scalar(n_dofs_ / n_components)
{
  // reinitialize
  reinitialize();
}

/**
 * \brief Reinitializes.
 */
template <int dim>
void DomainInvariantDiffusion<dim>::reinitialize()
{
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
void DomainInvariantDiffusion<dim>::compute_diffusion_matrix(
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

    // get normals and gradient norms
    auto gradient_norms = gradient_matrix->get_gradient_norms(i);
    auto normals = gradient_matrix->get_normals(i);
    auto column_indices = gradient_matrix->get_column_indices(i);

    // number of nonzero columns in row i
    const unsigned int n_columns = column_indices.size();

    // loop over row entries of scalar matrices
    for (unsigned int k = 0; k < n_columns; ++k)
    {
      // get column index
      const unsigned int j = column_indices[k];

      // compute indices of first component in vector system
      const unsigned int i0 = i * n_components;
      const unsigned int j0 = j * n_components;

      // make sure this is an off-diagonal entry
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
          // is no need to compute both max wave speeds here.

          // create tensor for normal vector
          Tensor<1, dim> normal_ji = gradient_matrix->get_normal(j, i);

          // compute maximum wave speed
          const double max_speed_ij =
            max_wave_speed->compute(solution_i, solution_j, normals[k]);
          const double max_speed_ji =
            max_wave_speed->compute(solution_j, solution_i, normal_ji);

          // get gradient norm
          const double c_ji = gradient_matrix->get_gradient_norm(j, i);

          // compute diffusion entry value (note D(i,j) = D(j,i))
          const double D_ij =
            -std::max(max_speed_ij * gradient_norms[k], max_speed_ji * c_ji);

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
    }   // end row iterator loop
  }     // end loop over rows
}
