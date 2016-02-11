/**
 * \file ShallowWaterBoundaryConditions.cc
 * \brief Provides the function definitions for the ShallowWaterBoundaryConditions class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] fe_ finite element system
 * \param[in] face_quadrature_ quadrature for face
 * \param[in] gravity_ acceleration due to gravity
 */
template <int dim>
ShallowWaterBoundaryConditions<dim>::ShallowWaterBoundaryConditions(
  const FESystem<dim> & fe_,
  const QGauss<dim - 1> & face_quadrature_,
  const double & gravity_)
  : BoundaryConditions<dim>(fe_, face_quadrature_),
    gravity(gravity_),
    height_extractor(0),
    momentum_extractor(1)
{
}

/**
 * \brief Adds face contributions to steady-state flux vector.
 *
 * Adds the following face contributions resulting from integration
 * by parts:
 * \f[
 *   r_i +=
 *     \left(\varphi_i^h,\mathbf{q}_h\cdot\mathbf{n}\right)_{\partial\Omega}
 *   + \left(\varphi_i^{\mathbf{q}},(\mathbf{q}\otimes\mathbf{v}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\cdot\mathbf{n}\right)_{\partial\Omega} .
 * \f]
 *
 * \param[in] height values at each quadrature point on face
 * \param[in] momentum values at each quadrature point on face
 * \param[inout] cell_residual cell steady-state residual, for which to
 *               add the face contributions
 */
template <int dim>
void ShallowWaterBoundaryConditions<dim>::integrate_face(
  const std::vector<double> & height,
  const std::vector<Tensor<1, dim>> & momentum,
  Vector<double> & cell_residual) const
{
  // identity tensor
  SymmetricTensor<2, dim> identity_tensor_sym = unit_symmetric_tensor<dim>();
  Tensor<2, dim> identity_tensor(identity_tensor_sym);

  // get normal vectors
  std::vector<Tensor<1, dim>> normal_vectors =
    this->fe_values_face.get_all_normal_vectors();

  for (unsigned int q = 0; q < this->n_quadrature_points_face; ++q)
  {
    // compute velocity
    Tensor<1, dim> velocity = momentum[q] / height[q];

    // compute inviscid flux for momentum
    Tensor<2, dim> momentum_inviscid_flux = outer_product(velocity, momentum[q]) +
      0.5 * gravity * std::pow(height[q], 2) * identity_tensor;

    // loop over DoFs in cell
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      cell_residual(i) +=
        (this->fe_values_face[height_extractor].value(i, q) * momentum[q] +
         this->fe_values_face[momentum_extractor].value(i, q) *
           momentum_inviscid_flux) *
        normal_vectors[q] * this->fe_values_face.JxW(q);
  }
}
