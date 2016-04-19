/**
 * \file ShallowWater2DDamBreakBC.cc
 * \brief Provides the function definitions for the ShallowWater2DDamBreakBC
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] fe_ finite element system
 * \param[in] face_quadrature_ quadrature for face
 * \param[in] gravity_ acceleration due to gravity
 */
template <int dim>
ShallowWater2DDamBreakBC<dim>::ShallowWater2DDamBreakBC(
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
void ShallowWater2DDamBreakBC<dim>::apply_boundary_condition(
  const Cell &,
  const FEValues<dim> &,
  const FEFaceValues<dim> & fe_values_face,
  const Vector<double> & solution,
  const double &,
  Vector<double> & cell_residual)
{
  // get solution values on face
  std::vector<double> height(this->n_quadrature_points_face);
  fe_values_face[this->height_extractor].get_function_values(solution, height);
  std::vector<Tensor<1, dim>> momentum(this->n_quadrature_points_face);
  fe_values_face[this->momentum_extractor].get_function_values(solution,
                                                               momentum);

  // call integrate face function
  this->integrate_face(height, momentum, cell_residual);
}

/**
 * \brief Applies interior boundary condition for a face.
 *
 * \param[in] cell cell iterator
 * \param[in] fe_values_cell FE values for cell
 * \param[in] fe_values_face FE values for face
 * \param[in] solution solution vector
 * \param[in] dt time step size \f$\Delta t\f$
 * \param[inout] cell_residual steady-state residual for cell
 */
template <int dim>
void ShallowWater2DDamBreakBC<dim>::apply_interior_boundary_condition(
  const Cell &,
  const FEValues<dim> &,
  const FEFaceValues<dim> & fe_values_face,
  const Vector<double> & solution,
  const double &,
  Vector<double> & cell_residual)
{
  // get quadrature points on face
  std::vector<Point<dim>> points(this->n_quadrature_points_face);
  points = fe_values_face.get_quadrature_points();

  // get normal vectors
  std::vector<Tensor<1, dim>> normal_vectors =
    this->fe_values_face.get_all_normal_vectors();

  // loop over quadrature points to determine if face is along interior wall
  bool face_on_interior_wall = false;
  for (unsigned int q = 0; q < this->n_quadrature_points_face; ++q)
  {
    if (std::abs(points[q][0]) < 1.0e-10)
    {
      // if the left half
      if (normal_vectors[q][0] < 0.0)
      {

      if (points[q][1] <= 560.0 || points[q][1] >= 840.0)
      {
        face_on_interior_wall = true;
        break;
      }
      }
    }
  }

  // if the face is on interior wall, apply wall boundary condition
  if (face_on_interior_wall)
  {
  // get solution values on face
  std::vector<double> height(this->n_quadrature_points_face);
  fe_values_face[this->height_extractor].get_function_values(solution, height);
  std::vector<Tensor<1, dim>> momentum(this->n_quadrature_points_face);
  fe_values_face[this->momentum_extractor].get_function_values(solution,
                                                               momentum);

  // identity tensor
  SymmetricTensor<2, dim> identity_tensor_sym = unit_symmetric_tensor<dim>();
  Tensor<2, dim> identity_tensor(identity_tensor_sym);

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

}
