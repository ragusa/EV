/**
 * \file EntropyViscosity.cc
 * \brief Provides the function definitions for the EntropyViscosity class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] parameters_ parameters
 * \param[in] entropy_ pointer to entropy
 * \param[in] cell_diameter_ map of cell iterators to cell diameters
 * \param[in] dof_handler_ degree of freedom handler
 */
template <int dim>
EntropyViscosity<dim>::EntropyViscosity(
  const RunParameters<dim> & parameters_,
  const std::shared_ptr<Entropy<dim>> & entropy_,
  const CellMap & cell_diameter_,
  const FESystem<dim> & fe_,
  const DoFHandler<dim> & dof_handler_,
  const QGauss<dim> & cell_quadrature_,
  const QGauss<dim - 1> & face_quadrature_)
  : Viscosity<dim>(dof_handler_),
    entropy(entropy_),
    residual_coefficient(parameters_.entropy_residual_coef),
    jump_coefficient(parameters_.entropy_jump_coef),
    smoothing_weight(parameters_.entropy_viscosity_smoothing_weight),
    cell_diameter(&cell_diameter_),
    fe(&fe_),
    cell_quadrature(&cell_quadrature_),
    face_quadrature(&face_quadrature_),
    n_q_points_cell(cell_quadrature_.size()),
    n_q_points_face(face_quadrature_.size()),
    faces_per_cell(GeometryInfo<dim>::faces_per_cell)
{
}

/**
 * \brief Computes viscosity values.
 *
 * \param[in] new_solution new solution
 * \param[in] old_solution old solution
 * \param[in] dt time step size
 * \param[in] n time index
 */
template <int dim>
void EntropyViscosity<dim>::update(const Vector<double> & new_solution,
                                   const Vector<double> & old_solution,
                                   const double & dt,
                                   const unsigned int &)
{
  // reinitialize entropy with new solution: reintializes group FE values
  // and normalization parameters
  entropy->reinitialize(new_solution);

  // compute entropy viscosity for each cell
  Cell cell = this->dof_handler->begin_active(), endc = this->dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // compute entropy residuals at each quadrature point in cell
    const auto entropy_residual =
      compute_entropy_residual(new_solution, old_solution, dt, cell);

    // compute max entropy gradient jump in cell
    const double max_entropy_jump = compute_max_entropy_jump(cell);

    // compute normalization at each quadrature point in cell
    const auto entropy_normalization =
      entropy->compute_entropy_normalization(new_solution, cell);

    // compute entropy viscosity
    const double h2 = std::pow(cell_diameter->at(cell), 2);
    this->values[cell] = 0.0;
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      this->values[cell] =
        std::max(this->values[cell],
                 h2 * (residual_coefficient * std::abs(entropy_residual[q]) +
                       jump_coefficient * max_entropy_jump) /
                   entropy_normalization[q]);
    }
  }
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
std::vector<double> EntropyViscosity<dim>::compute_entropy_residual(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const Cell & cell) const
{
  // FE values
  FEValues<dim> fe_values(
    *fe, *cell_quadrature, update_values | update_gradients);
  fe_values.reinit(cell);

  // compute entropy of current and old solutions
  auto entropy_new = entropy->compute_entropy(new_solution, fe_values);
  auto entropy_old = entropy->compute_entropy(old_solution, fe_values);
  auto divergence_entropy_flux = entropy->compute_divergence_entropy_flux(cell);

  // compute entropy residual at each quadrature point on cell
  std::vector<double> entropy_residual(n_q_points_cell);
  for (unsigned int q = 0; q < n_q_points_cell; ++q)
  {
    entropy_residual[q] =
      (entropy_new[q] - entropy_old[q]) / dt + divergence_entropy_flux[q];
  }

  return entropy_residual;
}

/**
 * \brief Computes the max entropy jump in a cell.
 *
 * \param[in] cell cell iterator
 *
 * \return max entropy jump in cell
 */
template <int dim>
double EntropyViscosity<dim>::compute_max_entropy_jump(const Cell & cell) const
{
  double max_jump_in_cell = 0.0;

  // loop over faces
  for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
  {
    // determine if face is interior
    typename DoFHandler<dim>::face_iterator face = cell->face(iface);
    if (face->at_boundary() == false)
    {
      // get gradients and normal vectors from this cell
      auto entropy_flux_gradients_this_cell =
        entropy->compute_entropy_flux_gradients_face(cell, iface);
      auto normal_vectors = entropy->get_normal_vectors(cell, iface);

      // get iterator to neighboring cell and face index of this face
      Assert(cell->neighbor(iface).state() == IteratorState::valid,
             ExcInternalError());
      Cell neighbor = cell->neighbor(iface);
      const unsigned int iface_neighbor = cell->neighbor_of_neighbor(iface);
      Assert(iface_neighbor < faces_per_cell, ExcInternalError());

      // get gradients from neighboring cell
      auto entropy_flux_gradients_neighbor_cell =
        entropy->compute_entropy_flux_gradients_face(neighbor, iface_neighbor);

      // loop over face quadrature points to determine max jump on face
      double max_jump_on_face = 0.0;
      for (unsigned int q = 0; q < n_q_points_face; ++q)
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

/**
 * \brief Smooths the entropy visosity profile using the maximum of the
 *        cell and its neighbors.
 */
template <int dim>
void EntropyViscosity<dim>::smooth_entropy_viscosity_max()
{
  // copy entropy viscosities
  CellMap entropy_viscosity_unsmoothed = this->values;

  // loop over cells
  Cell cell = this->dof_handler->begin_active(), endc = this->dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // loop over faces in cell
    for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
    {
      // determine if face is interior
      typename DoFHandler<dim>::face_iterator face = cell->face(iface);
      if (face->at_boundary() == false)
      {
        // get the cell's neighbor through this interior face
        Cell neighbor_cell = cell->neighbor(iface);

        // take max of cell and its neighbor
        this->values[cell] = std::max(
          this->values[cell], entropy_viscosity_unsmoothed[neighbor_cell]);
      }
    }
  }
}

/**
 * \brief Smooths the entropy visosity profile using a weighted average of the
 *        cell and its neighbors.
 */
template <int dim>
void EntropyViscosity<dim>::smooth_entropy_viscosity_average()
{
  // copy entropy viscosities
  CellMap entropy_viscosity_unsmoothed = this->values;

  // loop over cells
  Cell cell = this->dof_handler->begin_active(), endc = this->dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // initialize sum for average
    double sum = smoothing_weight * entropy_viscosity_unsmoothed[cell];

    // initialize neighbor count
    unsigned int neighbor_count = 0;

    // loop over faces in cell
    for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
    {
      // determine if face is interior
      typename DoFHandler<dim>::face_iterator face = cell->face(iface);
      if (face->at_boundary() == false)
      {
        // get the cell's neighbor through this interior face
        Cell neighbor_cell = cell->neighbor(iface);

        // add to sum
        sum += entropy_viscosity_unsmoothed[neighbor_cell];

        // increment neighbor count
        neighbor_count++;
      }
    }

    // compute average, smoothed value
    this->values[cell] = sum / (neighbor_count + smoothing_weight);
  }
}
