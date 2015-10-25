/**
 * \brief Constructor.
 *
 * \param[in] fe_ finite element system
 * \param[in] face_quadrature_ quadrature for face
 * \param[in] gravity_ acceleration due to gravity
 */
template <int dim>
ShallowWaterWallBC<dim>::ShallowWaterWallBC(
  const FESystem<dim> & fe_,
  const QGauss<dim - 1> & face_quadrature_,
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
 * \param[in] dt time step size \f$\Delta t\f$
 * \param[inout] cell_residual steady-state residual for cell
 */
template <int dim>
void ShallowWaterWallBC<dim>::apply_boundary_condition(
  const Cell &,
  const FEValues<dim> &,
  const FEFaceValues<dim> & fe_values_face,
  const Vector<double> & solution,
  const double &,
  Vector<double> & cell_residual)
{
  // get height values on face
  std::vector<double> height(this->n_quadrature_points_face);
  fe_values_face[this->height_extractor].get_function_values(solution, height);

  // use zero momentum values
  std::vector<Tensor<1, dim>> momentum(this->n_quadrature_points_face);

  // call integrate face function
  this->integrate_face(height, momentum, cell_residual);
}
