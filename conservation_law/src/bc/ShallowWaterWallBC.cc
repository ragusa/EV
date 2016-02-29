/**
 * \file ShallowWaterWallBC.cc
 * \brief Provides the function definitions for the ShallowWaterWallBC class.
 */

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

  // identity tensor
  SymmetricTensor<2, dim> identity_tensor_sym = unit_symmetric_tensor<dim>();
  Tensor<2, dim> identity_tensor(identity_tensor_sym);

  // get normal vectors
  std::vector<Tensor<1, dim>> normal_vectors =
    this->fe_values_face.get_all_normal_vectors();

  for (unsigned int q = 0; q < this->n_quadrature_points_face; ++q)
  {
    // compute inviscid flux for momentum with v.n = 0
    Tensor<2, dim> momentum_inviscid_flux =
      0.5 * this->gravity * std::pow(height[q], 2) * identity_tensor;

    // loop over DoFs in cell
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      cell_residual(i) +=
        this->fe_values_face[this->momentum_extractor].value(i, q) *
        momentum_inviscid_flux * normal_vectors[q] * this->fe_values_face.JxW(q);
  }
}
