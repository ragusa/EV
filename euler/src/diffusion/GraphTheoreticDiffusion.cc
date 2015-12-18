/**
 * \file GraphTheoreticDiffusion.cc
 * \brief Provides the function definitions for the GraphTheoreticDiffusion class.
 */

/**
 * \brief Constructor.
 *
 * Here, the unit cell matrix for local bilinear forms, \f$\mathrm{B}^K\f$,
 * is computed:
 * \f[
 *   \mathrm{B}^K_{i,j} = \frac{b_K(\varphi_i,\varphi_j)}{|K|} \,,
 * \f]
 * where \f$b_K(\varphi_i,\varphi_j)\f$ is the local viscous bilinear form:
 * \f[
 *   b_K(\varphi_i,\varphi_j) = \left\{\begin{array}{c l}
 *     |K|                & j = i, \quad i,j\in\mathcal{I}(K) \\
 *     -\frac{|K|}{n_K-1} & j \ne i, \quad i,j\in\mathcal{I}(K), \quad m(j)=m(i)
 * \\
 *     0 & \mbox{otherwise}
 *   \end{array}\right. \,,
 * \f]
 * where \f$|K|\f$ is the cell volume, \f$n_K=\mbox{card}(\mathcal{I}(K))\f$,
 * and \f$m(i)\f$ is the component index of degree of freedom \f$i\f$.
 *
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] dofs_per_cell_ number of degrees of freedom per cell
 * \param[in] n_components_ number of solution components
 */
template <int dim>
GraphTheoreticDiffusion<dim>::GraphTheoreticDiffusion(
  const DoFHandler<dim> & dof_handler_,
  const unsigned int & dofs_per_cell_,
  const unsigned int & n_components_)
  : ArtificialDiffusion<dim>(),
    dof_handler(&dof_handler_),
    dofs_per_cell(dofs_per_cell_),
    n_components(n_components_),
    dofs_per_cell_per_component(dofs_per_cell_ / n_components_),
    unit_cell_matrix(dofs_per_cell, dofs_per_cell)
{
  // compute unit cell matrix
  for (unsigned int m = 0; m < n_components; ++m)
    for (unsigned int i = 0; i < dofs_per_cell_per_component; ++i)
    {
      const unsigned int ii = i * n_components + m;
      for (unsigned int j = 0; j < dofs_per_cell_per_component; ++j)
      {
        const unsigned int jj = j * n_components + m;
        if (ii == jj)
          unit_cell_matrix(ii, jj) = 1.0;
        else
          unit_cell_matrix(ii, jj) = -1.0 / (dofs_per_cell_per_component - 1.0);
      }
    }
}

/**
 * \brief Computes artificial diffusion matrix.
 *
 * The artificial diffusion matrix is computed as
 * \f[
 *   D_{i,j} = \sum\limits_{K\subset S_{i,j}} \nu_K b_K(\varphi_i, \varphi_j) \,,
 * \f]
 * where \f$\nu_K\f$ is the cell viscosity.
 *
 * \param[in] solution solution vector
 * \param[in] viscosity pointer to viscosity cell map
 * \param[out] diffusion_matrix diffusion matrix
 */
template <int dim>
void GraphTheoreticDiffusion<dim>::compute_diffusion_matrix(
  const Vector<double> &,
  const std::shared_ptr<Viscosity<dim>> viscosity,
  SparseMatrix<double> & diffusion_matrix)
{
  // reset matrix
  diffusion_matrix = 0;

  // local DoF indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler->begin_active(),
                                                 endc = dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // get global degree of freedom indices
    cell->get_dof_indices(local_dof_indices);

    // get cell volume
    const double cell_volume = cell->measure();

    // aggregate local contribution into global matrix
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        diffusion_matrix.add(local_dof_indices[i],
                             local_dof_indices[j],
                             (*viscosity)[cell] * cell_volume *
                               unit_cell_matrix(i, j));
  }
}
