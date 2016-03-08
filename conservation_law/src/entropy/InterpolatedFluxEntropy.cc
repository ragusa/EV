/**
 * \file InterpolatedFluxEntropy.cc
 * \brief Provides the function definitions for the InterpolatedFluxEntropy class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] use_max_entropy_deviation_normalization_ flag to use max
 *            entropy deviation from average as entropy normalization
 * \param[in] domain_volume_ volume of domain
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] fe_ finite element system
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] face_quadrature_ face quadrature
 */
template <int dim>
InterpolatedFluxEntropy<dim>::InterpolatedFluxEntropy(
  const bool & use_max_entropy_deviation_normalization_,
  const double & domain_volume_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const QGauss<dim> & cell_quadrature_,
  const QGauss<dim - 1> & face_quadrature_)
  : Entropy<dim>(use_max_entropy_deviation_normalization_,
                 domain_volume_,
                 dof_handler_,
                 fe_,
                 cell_quadrature_,
                 face_quadrature_)
{
}

/**
 * \brief Computes the entropy residual at each quadrature point in a cell.
 *
 * \param[in] new_solution new solution
 * \param[in] old_solution old solution
 * \param[in] dt time step size
 * \param[in] cell cell iterator
 *
 * \return vector of entropy residual for each quadrature point in cell
 */
template <int dim>
std::vector<double> InterpolatedFluxEntropy<dim>::compute_entropy_residual(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const Cell & cell)
{
  // FE values
  FEValues<dim> fe_values(
    *this->fe, *this->cell_quadrature, update_values | update_gradients);
  fe_values.reinit(cell);

  // compute entropy of current and old solutions
  std::vector<double> entropy_new = compute_entropy(new_solution, fe_values);
  std::vector<double> entropy_old = compute_entropy(old_solution, fe_values);
  std::vector<double> divergence_entropy_flux =
    compute_divergence_entropy_flux(cell);

  // compute entropy residual at each quadrature point on cell
  std::vector<double> entropy_residual(this->n_q_points_cell);
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
  {
    entropy_residual[q] =
      (entropy_new[q] - entropy_old[q]) / dt + divergence_entropy_flux[q];
  }

  return entropy_residual;
}

/**
 * \brief Computes the max entropy jump in a cell.
 *
 * The max entropy jump on a cell is computed as
 * \f[
 *   J_K = \max\limits_{F\in\mathcal{F}(K)} J_F \,.
 * \f]
 *
 * \param[in] solution  solution vector
 * \param[in] cell cell iterator
 *
 * \return max entropy jump in cell
 */
template <int dim>
double InterpolatedFluxEntropy<dim>::compute_max_entropy_jump(
  const Vector<double> &, const Cell & cell)
{
  // initialize max entropy jump to zero
  double max_jump_in_cell = 0.0;

  // loop over faces
  for (unsigned int iface = 0; iface < this->faces_per_cell; ++iface)
  {
    // determine if face is interior
    typename DoFHandler<dim>::face_iterator face = cell->face(iface);
    if (face->at_boundary() == false)
    {
      // get gradients and normal vectors from this cell
      auto entropy_flux_gradients_this_cell =
        compute_entropy_flux_gradients_face(cell, iface);
      auto normal_vectors = get_normal_vectors(cell, iface);

      // get iterator to neighboring cell and face index of this face
      Assert(cell->neighbor(iface).state() == IteratorState::valid,
             ExcInternalError());
      Cell neighbor = cell->neighbor(iface);
      const unsigned int iface_neighbor = cell->neighbor_of_neighbor(iface);
      Assert(iface_neighbor < this->faces_per_cell, ExcInternalError());

      // get gradients from neighboring cell
      auto entropy_flux_gradients_neighbor_cell =
        compute_entropy_flux_gradients_face(neighbor, iface_neighbor);

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
