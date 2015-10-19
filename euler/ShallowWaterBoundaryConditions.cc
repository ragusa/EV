/**
 * \brief Constructor.
 */
template <int dim>
ShallowWaterBoundaryConditions<dim>::ShallowWaterBoundaryConditions()
{
}

/**
 * \brief
 */
template <int dim>
void ShallowWaterBoundaryConditions<dim>::integrate_face(
  const std::vector<double> & height,
  const std::vector<Tensor<1,dim>> & momentum
  Vector<double> & cell_residual) const
{
  for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      cell_residual(i) += (
        fe_values_face[height_extractor].value(i,q)
          * momentum[q] * normal_vectors[q]
        + fe_values_face[momentum_extractor]
          * momentum_inviscid_flux[q]
        )*fe_values_face.JxW(q);
  }
}
