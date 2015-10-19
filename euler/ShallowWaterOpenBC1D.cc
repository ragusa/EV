/**
 * \brief Constructor.
 *
 * \param[in] fe_ finite element system
 * \param[in] face_quadrature_ quadrature for face
 * \param[in] gravity_ acceleration due to gravity
 */
template <int dim>
ShallowWaterOpenBC1D<dim>::ShallowWaterOpenBC1D(
  const FESystem<dim> & fe_,
  const QGauss<dim> & face_quadrature_,
  const double & gravity_)
  : ShallowWaterBoundaryConditions<dim>(fe_, face_quadrature_, gravity_)
{
}

/**
 * \brief Applies boundary condition for a face.
 *
 * \param[in] cell cell iterator
 * \param[in] fe_values_cell FE values for cell
 * \param[in] fe_values_face FE values for face
 * \param[in] solution solution vector
 * \param[inout] cell_residual steady-state residual for cell
 */
void ShallowWaterOpenBC1D<dim>::apply_boundary_condition(
  const Cell & cell,
  const FEValues<dim> & fe_values_cell,
  const FEValuesFace<dim> & fe_values_face,
  const Vector<double> & solution,
  Vector<double> & cell_residual)
{
  // compute height and momentum values
  std::vector<double> height(n_quadrature_points_face, 0.0);
  std::vector<Tensor<1, dim>> momentum(n_quadrature_points_face, 0.0);

  // call integrate face function
  integrate_face(height, momentum, cell_residual);
}
