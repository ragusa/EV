/**
 * \file ScalarEntropy.cc
 * \brief Provides the function definitions for the ScalarEntropy class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] domain_volume_ volume of domain
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] fe_ finite element system
 * \param[in] triangulation_ triangulation
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] face_quadrature_ face quadrature
 */
template <int dim>
ScalarEntropy<dim>::ScalarEntropy(const double & domain_volume_,
                                  const DoFHandler<dim> & dof_handler_,
                                  const FESystem<dim> & fe_,
  const Triangulation<dim> & triangulation_,
                                  const QGauss<dim> & cell_quadrature_,
                                  const QGauss<dim - 1> & face_quadrature_)
  : Entropy<dim>(
      true, domain_volume_, dof_handler_, fe_, cell_quadrature_, face_quadrature_),
    entropy_flux_fe_values_cell(dof_handler_, triangulation_, cell_quadrature_),
    entropy_flux_fe_values_face(dof_handler_, triangulation_, face_quadrature_)
{
}

/**
 * \brief Computes entropy \f$\eta\f$ at each quadrature point on cell or face.
 *
 * For now, the entropy function is assumed to be \f$\eta(u) = \frac{1}{2}u^2\f$.
 *
 * \param[in] solution solution vector
 * \param[in] fe_values FE values, already initialized for cell
 *
 * \return vector of entropy values at each quadrature point
 */
template <int dim>
std::vector<double> ScalarEntropy<dim>::compute_entropy(
  const Vector<double> & solution, const FEValuesBase<dim> & fe_values) const
{
  // get number of quadrature points
  const unsigned int n = fe_values.n_quadrature_points;

  // get solution at quadrature points
  std::vector<double> solution_local(n);
  fe_values.get_function_values(solution, solution_local);

  // compute entropy
  std::vector<double> entropy(n);
  for (unsigned int q = 0; q < n; ++q)
    entropy[q] = 0.5 * solution_local[q] * solution_local[q];

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
std::vector<double> ScalarEntropy<dim>::compute_divergence_entropy_flux(
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
std::vector<Tensor<2, dim>> ScalarEntropy<
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
std::vector<Tensor<1, dim>> ScalarEntropy<dim>::get_normal_vectors(
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
std::vector<double> ScalarEntropy<dim>::compute_entropy_normalization(
  const Vector<double> & solution, const Cell & cell) const
{
  // call function pointer
  return this->compute_max_entropy_deviation_normalization(solution, cell);
}

/**
 * \brief Reinitializes the entropy flux FE values.
 *
 * \param[in] solution solution vector
 */
template <int dim>
void ScalarEntropy<dim>::reinitialize_group_fe_values(
  const Vector<double> & solution)
{
  entropy_flux_fe_values_cell.reinitialize(solution);
  entropy_flux_fe_values_face.reinitialize(solution);
}
