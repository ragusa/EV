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
  const ShallowWaterRunParameters & parameters_,
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
  : InterpolatedFluxEntropy<dim>(parameters_.entropy_normalization == "average",
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
    AssertThrow(false, ExcNotImplemented());
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
 * \brief Computes the local entropy viscosity normalization coefficient
 *        \f$gh^2\f$ for each quadrature point in a cell.
 *
 * The local entropy viscosity normalization coefficient is computed as
 *
 * \f[
 *   \hat{\eta}_q = g h_q^2
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
