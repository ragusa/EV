/**
 * \file LaplacianDiffusion.cc
 * \brief Provides the function definitions for the LaplacianDiffusion class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] scalar_extractors_ vector of scalar extractors
 * \param[in] vector_extractors_ vector of vector extractors
 * \param[in] viscosity_ pointer to viscosity
 * \param[in] n_quadrature_points_ number of quadrature points in cell
 * \param[in] dofs_per_cell_ number of degrees of freedom per cell
 */
template <int dim>
LaplacianDiffusion<dim>::LaplacianDiffusion(
  const std::vector<FEValuesExtractors::Scalar> & scalar_extractors_,
  const std::vector<FEValuesExtractors::Vector> & vector_extractors_,
  const unsigned int & n_quadrature_points_,
  const unsigned int & dofs_per_cell_)
  : ArtificialDiffusion<dim>(),
    scalar_extractors(scalar_extractors_),
    vector_extractors(vector_extractors_),
    n_scalar_components(scalar_extractors_.size()),
    n_vector_components(vector_extractors_.size()),
    n_quadrature_points(n_quadrature_points_),
    dofs_per_cell(dofs_per_cell_)
{
}

/**
 * \brief Applies Laplacian diffusion for a cell to its residual.
 *
 * \param[in] viscosity viscosity
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 * \param[in] fe_values FE values, reinitialized for cell already
 * \param[inout] cell_residual the residual vector for the cell
 */
template <int dim>
void LaplacianDiffusion<dim>::apply(std::shared_ptr<Viscosity<dim>> viscosity,
                                    const Vector<double> & solution,
                                    const Cell & cell,
                                    const FEValues<dim> & fe_values,
                                    Vector<double> & cell_residual) const
{
  // get solution gradients for scalar components
  std::vector<std::vector<Tensor<1, dim>>> scalar_gradients(
    n_scalar_components, std::vector<Tensor<1, dim>>(n_quadrature_points));
  for (unsigned int j = 0; j < n_scalar_components; ++j)
    fe_values[scalar_extractors[j]].get_function_gradients(solution,
                                                           scalar_gradients[j]);

  // get solution gradients for vector components
  std::vector<std::vector<Tensor<2, dim>>> vector_gradients(
    n_vector_components, std::vector<Tensor<2, dim>>(n_quadrature_points));
  for (unsigned int j = 0; j < n_vector_components; ++j)
    fe_values[vector_extractors[j]].get_function_gradients(solution,
                                                           vector_gradients[j]);

  // compute contributions to cell residual
  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      // initialize sum for DoF i, quadrature point q, to zero
      double sum_i_q = 0.0;

      // loop over scalar components
      for (unsigned int j = 0; j < n_scalar_components; ++j)
        sum_i_q += fe_values[scalar_extractors[j]].gradient(i, q) *
          (-(*viscosity)[cell] * scalar_gradients[j][q]);

      // loop over vector components
      for (unsigned int j = 0; j < n_vector_components; ++j)
        sum_i_q += double_contract<0, 0, 1, 1>(
          fe_values[vector_extractors[j]].gradient(i, q),
          -(*viscosity)[cell] * vector_gradients[j][q]);

      // multiply by Jacobian
      sum_i_q *= fe_values.JxW(q);

      // add sum to cell residual
      cell_residual(i) += sum_i_q;
    }
  }
}
