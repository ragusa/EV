/**
 * \file ShallowWaterEntropy.cc
 * \brief Provides the function definitions for the ShallowWaterEntropy class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] parameters_ parameters
 * \param[in] height_extractor_ FE value extractor for height
 * \param[in] momentum_extractor_ FE value extractor for momentum
 * \param[in] gravity_ acceleration due to gravity
 * \param[in] bathymetry_vector_ bathymetry vector
 * \param[in] entropy_normalization_option_ entropy normalization option
 * \param[in] domain_volume_ volume of domain
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] fe_ finite element system
 * \param[in] triangulation_ triangulation
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] face_quadrature_ face quadrature
 */
template <int dim>
ShallowWaterEntropy<dim>::ShallowWaterEntropy(
  const ShallowWaterParameters<dim> & parameters_,
  const FEValuesExtractors::Scalar & height_extractor_,
  const FEValuesExtractors::Vector & momentum_extractor_,
  const double & gravity_,
  const Vector<double> & bathymetry_vector_,
  const double & domain_volume_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const Triangulation<dim> & triangulation_,
  const QGauss<dim> & cell_quadrature_,
  const QGauss<dim - 1> & face_quadrature_)
  : Entropy<dim>(parameters_.entropy_normalization == "average",
                 domain_volume_,
                 dof_handler_,
                 fe_,
                 cell_quadrature_,
                 face_quadrature_),
    compute_entropy_normalization_ptr(nullptr),
    height_extractor(height_extractor_),
    momentum_extractor(momentum_extractor_),
    gravity(gravity_),
    entropy_flux_fe_values_cell(dof_handler_,
                                triangulation_,
                                cell_quadrature_,
                                bathymetry_vector_,
                                gravity_),
    entropy_flux_fe_values_face(dof_handler_,
                                triangulation_,
                                face_quadrature_,
                                bathymetry_vector_,
                                gravity_),
    entropy_normalization_option(parameters_.entropy_normalization)
{
  // determine which function to use to compute entropy normalization
  if (entropy_normalization_option == "average")
    compute_entropy_normalization_ptr =
      &ShallowWaterEntropy<dim>::compute_max_entropy_deviation_normalization;
  else if (entropy_normalization_option == "local")
    compute_entropy_normalization_ptr =
      &ShallowWaterEntropy<dim>::compute_local_entropy_normalization;
  else
    ExcNotImplemented();
}

/**
 * \brief Computes entropy \f$\eta\f$ at each quadrature point on cell or face.
 *
 * For the shallow water equations, the entropy is defined as
 * \f[
 *   \eta(\mathbf{u}) = \frac{1}{2}\frac{\mathbf{q}\cdot\mathbf{q}}{h}
 *   + \frac{1}{2}g h^2
 * \f]
 */
template <int dim>
std::vector<double> ShallowWaterEntropy<dim>::compute_entropy(
  const Vector<double> & solution, const FEValuesBase<dim> & fe_values) const
{
  // get number of quadrature points
  const unsigned int n = fe_values.n_quadrature_points;

  // get height and momentum
  std::vector<double> height(n);
  std::vector<Tensor<1, dim>> momentum(n);
  fe_values[height_extractor].get_function_values(solution, height);
  fe_values[momentum_extractor].get_function_values(solution, momentum);

  // compute entropy
  std::vector<double> entropy(n);
  for (unsigned int q = 0; q < n; ++q)
    entropy[q] = 0.5 * momentum[q] * momentum[q] / height[q] +
      0.5 * gravity * height[q] * height[q];

  return entropy;
}

/**
 * \brief Computes divergence of entropy flux at each quadrature point in cell.
 *
 * This function uses group FEM values for entropy flux to interpolate values.
 *
 * \param[in] cell cell iterator
 *
 * \return vector of divergence of entropy flux at each quadrature point in cell.
 */
template <int dim>
std::vector<double> ShallowWaterEntropy<dim>::compute_divergence_entropy_flux(
  const Cell & cell)
{
  // reinitialize values for cell
  entropy_flux_fe_values_cell.reinit(cell);

  // get divergence of entropy flux
  auto divergence_entropy_flux =
    entropy_flux_fe_values_cell.get_function_divergences();

  return divergence_entropy_flux;
}

/**
 * \brief Computes gradients of entropy flux at each quadrature point in face.
 *
 * This function uses group FEM values for entropy flux to interpolate values.
 *
 * \param[in] cell cell iterator
 * \param[in] i_face index of face on cell
 *
 * \return vector of gradients of entropy flux at each quadrature point in cell.
 */
template <int dim>
std::vector<Tensor<2, dim>> ShallowWaterEntropy<
  dim>::compute_entropy_flux_gradients_face(const Cell & cell,
                                            const unsigned int & i_face)
{
  // reinitialize values for face
  entropy_flux_fe_values_face.reinit(cell, i_face);

  // get gradients of entropy flux
  std::vector<Tensor<2, dim>> entropy_flux_gradients(this->n_q_points_face);
  entropy_flux_fe_values_face.get_function_gradients(entropy_flux_gradients);

  return entropy_flux_gradients;
}

/**
 * \brief Gets normal vectors at each quadrature point in face.
 *
 * \param[in] cell cell iterator
 * \param[in] i_face index of face on cell
 *
 * \return vector of normal vectors at each quadrature point in face
 */
template <int dim>
std::vector<Tensor<1, dim>> ShallowWaterEntropy<dim>::get_normal_vectors(
  const Cell & cell, const unsigned int & i_face)
{
  // reinitialize values for face
  entropy_flux_fe_values_face.reinit(cell, i_face);

  return entropy_flux_fe_values_face.get_normal_vectors();
}

/**
 * \brief Computes the entropy viscosity normalization coefficient
 *        for each quadrature point in a cell using a function pointer.
 *
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 *
 * \return vector of entropy viscosity normalization coefficient for each
 *         quadrature point in cell
 */
template <int dim>
std::vector<double> ShallowWaterEntropy<dim>::compute_entropy_normalization(
  const Vector<double> & solution, const Cell & cell) const
{
  // call function pointer
  return (this->*compute_entropy_normalization_ptr)(solution, cell);
}

/**
 * \brief Computes the entropy residual at each quadrature point in a cell.
 *
 * \param[in] new_solution new solution
 * \param[in] old_solution old solution
 * \param[in] entropy_flux_fe_values FE values for entropy flux
 * \param[in] dt time step size
 * \param[in] cell cell iterator
 *
 * \return vector of entropy residual for each quadrature point in cell
 */
/*
template <int dim>
std::vector<double> ShallowWater<dim>::compute_entropy_residual(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const ShallowWaterEntropyFluxFEValuesCell<dim> & entropy_flux_fe_values,
  const double & dt,
  const Cell & cell) const
{
  // FE values
  FEValues<dim> fe_values(
    this->fe, this->cell_quadrature, update_values | update_gradients);
  fe_values.reinit(cell);

  Vector<double> entropy_new(this->n_q_points_cell);
  Vector<double> entropy_old(this->n_q_points_cell);
  std::vector<double> entropy_residual(this->n_q_points_cell);

  // compute entropy of current and old solutions
  compute_entropy(new_solution, fe_values, entropy_new);
  compute_entropy(old_solution, fe_values, entropy_old);
  std::vector<double> divergence_entropy_flux =
    entropy_flux_fe_values.get_function_divergences();

  // compute entropy residual at each quadrature point on cell
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
  {
    // compute entropy residual
    double dsdt = (entropy_new[q] - entropy_old[q]) / dt;
    entropy_residual[q] = dsdt + divergence_entropy_flux[q];
  }

  return entropy_residual;
}
*/

/**
 * \brief Computes the max entropy jump in a cell.
 *
 * \param[in] fe_values entropy flux FE values
 * \param[in] cell cell iterator
 *
 * \return max entropy jump in cell
 */
/*
template <int dim>
double ShallowWater<dim>::compute_max_entropy_jump(
  ShallowWaterEntropyFluxFEValuesFace<dim> & fe_values, const Cell & cell) const
{
  std::vector<Tensor<2, dim>> entropy_flux_gradients_this_cell(
    this->n_q_points_face);
  std::vector<Tensor<2, dim>> entropy_flux_gradients_neighbor_cell(
    this->n_q_points_face);
  std::vector<Tensor<1, dim>> normal_vectors(this->n_q_points_face);

  double max_jump_in_cell = 0.0;

  // loop over faces
  for (unsigned int iface = 0; iface < this->faces_per_cell; ++iface)
  {
    // determine if face is interior
    typename DoFHandler<dim>::face_iterator face = cell->face(iface);
    if (face->at_boundary() == false)
    {
      // reinitialize FE values
      fe_values.reinit(cell, iface);

      // get gradients
      fe_values.get_function_gradients(entropy_flux_gradients_this_cell);

      // get normal vectors
      normal_vectors = fe_values.get_normal_vectors();

      // get iterator to neighboring cell and face index of this face
      Assert(cell->neighbor(iface).state() == IteratorState::valid,
             ExcInternalError());
      Cell neighbor = cell->neighbor(iface);
      const unsigned int ineighbor = cell->neighbor_of_neighbor(iface);
      Assert(ineighbor < this->faces_per_cell, ExcInternalError());

      // get gradients from neighboring cell
      fe_values.reinit(neighbor, ineighbor);
      fe_values.get_function_gradients(entropy_flux_gradients_neighbor_cell);

      // loop over face quadrature points to determine max jump on face
      double max_jump_on_face = 0.0;
      for (unsigned int q = 0; q < this->n_q_points_face; ++q)
      {
        // compute difference in gradients across face
        Tensor<2, dim> entropy_flux_gradient_jump =
          entropy_flux_gradients_this_cell[q] -
          entropy_flux_gradients_neighbor_cell[q];
        double jump_on_face = std::abs(entropy_flux_gradient_jump *
                                       normal_vectors[q] * normal_vectors[q]);
        max_jump_on_face = std::max(max_jump_on_face, jump_on_face);
      }

      // update max jump in cell
      max_jump_in_cell = std::max(max_jump_in_cell, max_jump_on_face);
    }
  }

  return max_jump_in_cell;
}
*/

/**
 * \brief Computes the local entropy viscosity normalization coefficient
 *        \f$gh^2\f$ for each quadrature point in a cell.
 *
 * The local entropy viscosity normalization coefficient is computed as
 *
 * \f[
 *   c^{\mbox{normalization}}_q = g h_q^2
 * \f]
 *
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 *
 * \return vector of entropy viscosity normalization coefficient for each
 *         quadrature point in cell
 */
template <int dim>
std::vector<double> ShallowWaterEntropy<dim>::compute_local_entropy_normalization(
  const Vector<double> & solution, const Cell & cell) const
{
  // get height values
  FEValues<dim> fe_values(*(this->fe), *(this->cell_quadrature), update_values);
  fe_values.reinit(cell);
  std::vector<double> height(this->n_q_points_cell);
  fe_values[height_extractor].get_function_values(solution, height);

  // compute normalization at each quadrature point
  std::vector<double> normalization(this->n_q_points_cell);
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    normalization[q] = gravity * std::pow(height[q], 2);

  return normalization;
}

/**
 * \brief Reinitializes the entropy flux FE values.
 *
 * \param[in] solution solution vector
 */
template <int dim>
void ShallowWaterEntropy<dim>::reinitialize_group_fe_values(
  const Vector<double> & solution)
{
  entropy_flux_fe_values_cell.reinitialize(solution);
  entropy_flux_fe_values_face.reinitialize(solution);
}
