/**
 * \file LaplacianDiffusion.cc
 * \brief Provides the function definitions for the LaplacianDiffusion class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] scalar_extractors_ vector of scalar extractors
 * \param[in] vector_extractors_ vector of vector extractors
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] fe_ finite element system
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] dofs_per_cell_ number of degrees of freedom per cell
 */
template <int dim>
LaplacianDiffusion<dim>::LaplacianDiffusion(
  const std::vector<FEValuesExtractors::Scalar> & scalar_extractors_,
  const std::vector<FEValuesExtractors::Vector> & vector_extractors_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const QGauss<dim> & cell_quadrature_,
  const unsigned int & dofs_per_cell_)
  : ArtificialDiffusion<dim>(),
    scalar_extractors(scalar_extractors_),
    vector_extractors(vector_extractors_),
    dof_handler(&dof_handler_),
    fe(&fe_),
    cell_quadrature(&cell_quadrature_),
    n_scalar_components(scalar_extractors_.size()),
    n_vector_components(vector_extractors_.size()),
    n_quadrature_points(cell_quadrature_.size()),
    dofs_per_cell(dofs_per_cell_)
{
}

/**
 * \brief Computes artificial diffusion matrix.
 *
 * The artificial diffusion matrix is computed as
 * \f[
 *   D_{i,j} = -\sum\limits_{K\subset S_{i,j}}
 *     \int\limits_K\varphi_i(\mathbf{x})\nabla\cdot(\nu_K\nabla\varphi_j(\mathbf{x}))
 * dV \,,
 * \f]
 * which after integration by parts and dropping the boundary term, is
 * \f[
 *   D_{i,j} = \sum\limits_{K\subset S_{i,j}}
 *     \nu_K\int\limits_K\nabla\varphi_i(\mathbf{x})\cdot\nabla\varphi_j(\mathbf{x})
 * dV \,,
 * \f]
 * where \f$\nu_K\f$ is the cell viscosity.
 *
 * \param[in] solution solution vector
 * \param[in] viscosity pointer to viscosity cell map
 * \param[out] diffusion_matrix diffusion matrix
 */
template <int dim>
void LaplacianDiffusion<dim>::compute_diffusion_matrix(
  const Vector<double> &,
  const std::shared_ptr<Viscosity<dim>> viscosity,
  SparseMatrix<double> & diffusion_matrix)
{
  // reset matrix
  diffusion_matrix = 0;

  // local DoF indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // cell matrix
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  // FE values
  FEValues<dim> fe_values(
    *fe, *cell_quadrature, update_gradients | update_JxW_values);

  // loop over cells
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler->begin_active(),
                                                 endc = dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // reset cell matrix
    cell_matrix = 0;

    // reinitialize fe values for cell
    fe_values.reinit(cell);

    // get global degree of freedom indices
    cell->get_dof_indices(local_dof_indices);

    // compute cell matrix
    for (unsigned int q = 0; q < n_quadrature_points; ++q)
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          // initialize sum for this quadrature point
          double sum_i_q = 0.0;

          // loop over scalar components
          for (unsigned int k = 0; k < n_scalar_components; ++k)
            sum_i_q += fe_values[scalar_extractors[k]].gradient(i, q) *
              (*viscosity)[cell] * fe_values[scalar_extractors[k]].gradient(j, q);

          // loop over vector components
          for (unsigned int k = 0; k < n_vector_components; ++k)
            sum_i_q += double_contract<0, 0, 1, 1>(
              fe_values[vector_extractors[k]].gradient(i, q),
              (*viscosity)[cell] *
                fe_values[vector_extractors[k]].gradient(j, q));

          // finish cell matrix entry
          cell_matrix(i, j) += sum_i_q * fe_values.JxW(q);
        }
      }
    }

    // aggregate local contribution into global matrix
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        diffusion_matrix.add(
          local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
  }
}
