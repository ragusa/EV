/**
 * \brief Constructor.
 *
 * \param[in] fe_ finite element system
 * \param[in] face_quadrature_ quadrature for face
 * \param[in] gravity_ acceleration due to gravity
 */
template <int dim>
ShallowWaterSubcriticalWallBC1D<dim>::ShallowWaterSubcriticalWallBC1D(
  const FESystem<dim> & fe_,
  const QGauss<dim - 1> & face_quadrature_,
  const double & gravity_)
  : ShallowWaterBoundaryConditions<dim>(fe_, face_quadrature_, gravity_)
{
  // assert that this is 1-D
  Assert(dim == 1, ExcImpossibleInDim(dim));

  // assert that number of face quadrature points is 1 (as it should be in 1-D)
  Assert(this->n_quadrature_points_face == 1, ExcInvalidState());
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
template <int dim>
void ShallowWaterSubcriticalWallBC1D<dim>::apply_boundary_condition(
  const Cell & cell,
  const FEValues<dim> &,
  const FEFaceValues<dim> & fe_values_face,
  const Vector<double> & solution,
  Vector<double> & cell_residual)
{
  // get solution values on face
  std::vector<double> height(this->n_quadrature_points_face);
  fe_values_face[this->height_extractor].get_function_values(solution, height);
  std::vector<Tensor<1, dim>> momentum(this->n_quadrature_points_face);
  fe_values_face[this->momentum_extractor].get_function_values(solution,
                                                               momentum);

  // assert that flow is subcritical
  const Tensor<1, dim> velocity(momentum[0] / height[0]);
  const double speed = velocity.norm();
  const double speed_of_sound = std::sqrt(this->gravity * height[0]);
  const double froude_number = speed / speed_of_sound;
  Assert(froude_number < 1, ExcFlowNotSubcritical(froude_number));

  // get normal vectors
  std::vector<Tensor<1, dim>> normal_vectors =
    fe_values_face.get_all_normal_vectors();

  // get local degree of freedom indices
  std::vector<unsigned int> dof_indices(this->dofs_per_cell);
  cell->get_dof_indices(dof_indices);

  // set velocity boundary value from wall condition: u.n = 0
  const double velocity_bc = 0.0;

  double speed_of_sound_bc;

  // determine whether the boundary is the left or right
  if (normal_vectors[0][0] < 0.0) // left boundary
  {
    // get adjacent degree of freedom values from solution vector
    const double height_adjacent = solution[dof_indices[2]];
    const double momentum_adjacent = solution[dof_indices[3]];

    // compute adjacent velocity and speed of sound
    const double velocity_adjacent = momentum_adjacent / height_adjacent;
    const double speed_of_sound_adjacent =
      std::sqrt(this->gravity * height_adjacent);

    // compute boundary speed of sound
    speed_of_sound_bc =
      speed_of_sound_adjacent + 0.5 * (velocity_bc - velocity_adjacent);
  }
  else // right boundary
  {
    // get adjacent degree of freedom values from solution vector
    const double height_adjacent = solution[dof_indices[0]];
    const double momentum_adjacent = solution[dof_indices[1]];

    // compute adjacent velocity and speed of sound
    const double velocity_adjacent = momentum_adjacent / height_adjacent;
    const double speed_of_sound_adjacent =
      std::sqrt(this->gravity * height_adjacent);

    // compute boundary speed of sound
    speed_of_sound_bc =
      speed_of_sound_adjacent + 0.5 * (velocity_adjacent - velocity_bc);
  }

  // compute boundary height value
  const double height_bc = std::pow(speed_of_sound_bc, 2) / this->gravity;

  // compute boundary momentum
  Tensor<1, dim> momentum_bc;
  momentum_bc[0] = velocity_bc * height_bc;

  // convert values to vectors
  std::vector<double> height_vector(1, height_bc);
  std::vector<Tensor<1, dim>> momentum_vector(1, momentum_bc);

  // call integrate face function
  this->integrate_face(height_vector, momentum_vector, cell_residual);
}
