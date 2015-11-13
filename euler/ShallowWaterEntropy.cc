/**
 * \file ShallowWaterEntropy.cc
 * \brief Provides the function definitions for the ShallowWaterEntropy class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ShallowWaterEntropy<dim>::ShallowWaterEntropy()
  : Entropy<dim>()
{
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
void ShallowWater<dim>::compute_entropy(const Vector<double> & solution,
                                        const FEValuesBase<dim> & fe_values,
                                        Vector<double> & entropy) const
{
  // get number of quadrature points
  const unsigned int n = entropy.size();

  // get height and momentum
  std::vector<double> height(n);
  std::vector<Tensor<1, dim>> momentum(n);
  fe_values[height_extractor].get_function_values(solution, height);
  fe_values[momentum_extractor].get_function_values(solution, momentum);

  // compute entropy
  for (unsigned int q = 0; q < n; ++q)
    entropy(q) = 0.5 * momentum[q] * momentum[q] / height[q] +
      0.5 * gravity * height[q] * height[q];
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

/**
 * \brief Computes the max entropy jump in a cell.
 *
 * \param[in] fe_values entropy flux FE values
 * \param[in] cell cell iterator
 *
 * \return max entropy jump in cell
 */
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
