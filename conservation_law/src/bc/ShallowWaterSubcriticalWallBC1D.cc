/**
 * \file ShallowWaterSubcriticalWallBC1D.cc
 * \brief Provides the function definitions for the
 * ShallowWaterSubcriticalWallBC1D class.
 */

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
 * \param[in] dt time step size \f$\Delta t\f$
 * \param[inout] cell_residual steady-state residual for cell
 */
template <int dim>
void ShallowWaterSubcriticalWallBC1D<dim>::apply_boundary_condition(
  const Cell & cell,
  const FEValues<dim> &,
  const FEFaceValues<dim> & fe_values_face,
  const Vector<double> & solution,
  const double & dt,
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

  // get position of face
  std::vector<Point<dim>> face_points = fe_values_face.get_quadrature_points();
  const double face_position = face_points[0][0];

  double speed_of_sound_bc;

  // determine whether the boundary is the left or right
  if (normal_vectors[0][0] < 0.0) // left boundary
  {
    /*
        // get interior degree of freedom values from solution vector
        const double height_interior = solution[dof_indices[2]];
        const double momentum_interior = solution[dof_indices[3]];

        // compute interior velocity and speed of sound
        const double velocity_interior = momentum_interior / height_interior;
        const double speed_of_sound_interior =
          std::sqrt(this->gravity * height_interior);
    */

    // estimate exiting wave velocity
    const double wave_velocity = velocity[0] - speed_of_sound;

    // get the interior velocity and speed of sound values
    double speed_of_sound_interior, velocity_interior;
    get_interior_values(face_position,
                        wave_velocity,
                        dt,
                        cell,
                        solution,
                        speed_of_sound_interior,
                        velocity_interior);

    // compute boundary speed of sound
    speed_of_sound_bc =
      speed_of_sound_interior + 0.5 * (velocity_bc - velocity_interior);
  }
  else // right boundary
  {
    /*
        // get interior degree of freedom values from solution vector
        const double height_interior = solution[dof_indices[0]];
        const double momentum_interior = solution[dof_indices[1]];

        // compute interior velocity and speed of sound
        const double velocity_interior = momentum_interior / height_interior;
        const double speed_of_sound_interior =
          std::sqrt(this->gravity * height_interior);
    */

    // estimate exiting wave velocity
    const double wave_velocity = velocity[0] + speed_of_sound;

    // get the interior velocity and speed of sound values
    double speed_of_sound_interior, velocity_interior;
    get_interior_values(face_position,
                        wave_velocity,
                        dt,
                        cell,
                        solution,
                        speed_of_sound_interior,
                        velocity_interior);

    // compute boundary speed of sound
    speed_of_sound_bc =
      speed_of_sound_interior + 0.5 * (velocity_interior - velocity_bc);
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

/**
 * \brief Gets the solution values at the interior point for the exiting
 *        characteristic.
 *
 * \param[in] face_position x-position of boundary face
 * \param[in] wave_velocity velocity of exiting characteristic
 * \param[in] dt time step size \f$\Delta t\f$
 * \param[in] cell iterator for boundary cell
 * \param[in] solution solution vector
 * \param[out] speed_of_sound_interior speed of sound at the beginning position
 *             of the exiting characteristic
 * \param[out] velocity_interior velocity at the beginning position
 *             of the exiting characteristic
 */
template <int dim>
void ShallowWaterSubcriticalWallBC1D<dim>::get_interior_values(
  const double & face_position,
  const double & wave_velocity,
  const double & dt,
  const Cell & cell,
  const Vector<double> & solution,
  double & speed_of_sound_interior,
  double & velocity_interior) const
{
  // compute the position of the beginning of the characteristic to
  // be integrated
  const double interior_position = face_position + wave_velocity * dt;
  const double dx = cell->measure();
  double left_cell_position;
  if (wave_velocity < 0.0) // cell is at left boundary
    left_cell_position = face_position;
  else
    left_cell_position = face_position - dx;
  const double unit_cell_position = (interior_position - left_cell_position) / dx;
  const std::vector<Point<1>> interior_point(1, Point<1>(unit_cell_position));

  // create dummy quadrature for evaluating solution at characteristic
  // origin point
  Quadrature<dim> interior_quadrature(interior_point);

  // create FE values for evaluating solution at interior point
  FEValues<dim> fe_values_interior(this->fe, interior_quadrature, update_values);
  fe_values_interior.reinit(cell);

  // get solution at interior point
  std::vector<double> height_interior(1);
  std::vector<Tensor<1, dim>> momentum_interior(1);
  fe_values_interior[this->height_extractor].get_function_values(solution,
                                                                 height_interior);
  fe_values_interior[this->momentum_extractor].get_function_values(
    solution, momentum_interior);

  // compute adjacent velocity and speed of sound
  velocity_interior = momentum_interior[0][0] / height_interior[0];
  speed_of_sound_interior = std::sqrt(this->gravity * height_interior[0]);
}
