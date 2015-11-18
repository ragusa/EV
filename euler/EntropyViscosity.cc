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
  const ConservationLawParameters<dim> & parameters_,
  const std::shared_ptr<Entropy<dim>> & entropy_,
  const CellMap & cell_diameter_,
  const FESystem<dim> & fe_,
  const DoFHandler<dim> & dof_handler_)
  : Viscosity<dim>(dof_handler_),
    entropy(entropy_),
    residual_coefficient(parameters.entropy_residual_coef),
    jump_coefficient(parameters.entropy_jump_coef),
    cell_diameter(&cell_diameter_),
    n_q_points_cell(parameters_.n_quadrature_points),
    fe(&fe_)
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
    const std::vector<double> entropy_residual =
      compute_entropy_residual(new_solution, old_solution, dt, cell);

    // compute max entropy gradient jump in cell
    const double max_entropy_jump = compute_max_entropy_jump(new_solution, cell);

    // compute normalization at each quadrature point in cell
    const std::vector<double> entropy_normalization =
      compute_entropy_normalization(new_solution, entropy_average, cell);

    // compute entropy viscosity
    double h2 = std::pow((*cell_diameter)[cell], 2);
    this->values[cell] = 0.0;
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      this->values[cell] =
        std::max(this->values[cell],
                 h2 * (residual_coefficient * std::abs(entropy_residual[q]) +
                       jump_coefficient * max_entropy_jump) /
                   entropy_normalization[q]);
  }
}

/**
 * Computes the average of the entropy over the domain.
 *
 * \param[in] solution solution with which to calculate entropy for the
 *            normalization constant
 *
 * \return average entropy in domain
 */
template <int dim>
double ConservationLaw<dim>::compute_average_entropy(
  const Vector<double> & solution) const
{
  FEValues<dim> fe_values(
    *fe, cell_quadrature, update_values | update_JxW_values);

  Vector<double> entropy(n_q_points_cell);

  double domain_integral_entropy = 0.0;

  // loop over cells
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // compute entropy
    compute_entropy(solution, fe_values, entropy);

    // loop over quadrature points
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      // add contribution of quadrature point to entropy integral
      domain_integral_entropy += entropy[q] * fe_values.JxW(q);
    }
  }
  // domain-averaged entropy_values
  const double domain_averaged_entropy = domain_integral_entropy / domain_volume;

  return domain_averaged_entropy;
}

/**
 * \brief Computes the entropy deviation from the average at each quadrature
 *        point in a cell.
 *
 * \param[in] solution solution with which to calculate entropy for the
 *            normalization constant
 * \param[in] entropy_average average entropy in domain
 *
 * \return vector of entropy deviations at each quadrature point in cell
 */
template <int dim>
std::vector<double> ConservationLaw<dim>::compute_entropy_normalization(
  const Vector<double> & solution,
  const double & entropy_average,
  const Cell & cell) const
{
  FEValues<dim> fe_values(fe, cell_quadrature, update_values);

  Vector<double> entropy(n_q_points_cell);

  std::vector<double> normalization_constant(n_q_points_cell);

  // reinitialize FE values
  fe_values.reinit(cell);

  // compute entropy
  compute_entropy(solution, fe_values, entropy);

  // loop over quadrature points
  for (unsigned int q = 0; q < n_q_points_cell; ++q)
    normalization_constant[q] = std::abs(entropy[q] - entropy_average);

  return normalization_constant;
}

/**
 * \brief Computes the entropy residual at each quadrature point in a cell.
 *
 * \param[in] new_solution new solution
 * \param[in] old_solution old solution
 * \param[in] dt time step size
 * \param[in] cell cell iterator
 *
 * \return vector of entropy residuals at each quadrature point in a cell
 */
template <int dim>
std::vector<double> EntropyViscosity<dim>::compute_entropy_residual(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const Cell & cell) const
{
  // FE values
  FEValues<dim> fe_values(fe, cell_quadrature, update_values | update_gradients);
  fe_values.reinit(cell);

  Vector<double> entropy_new(n_q_points_cell);
  Vector<double> entropy_old(n_q_points_cell);
  Vector<double> divergence_entropy_flux(n_q_points_cell);
  std::vector<double> entropy_residual(n_q_points_cell);

  // compute entropy of current and old solutions
  compute_entropy(new_solution, fe_values, entropy_new);
  compute_entropy(old_solution, fe_values, entropy_old);
  compute_divergence_entropy_flux(
    new_solution, fe_values, divergence_entropy_flux);

  // compute entropy residual at each quadrature point on cell
  for (unsigned int q = 0; q < n_q_points_cell; ++q)
  {
    // compute entropy residual
    double dsdt = (entropy_new[q] - entropy_old[q]) / dt;
    entropy_residual[q] = std::abs(dsdt + divergence_entropy_flux[q]);
  }

  return entropy_residual;
}

/**
 * \brief Computes the max entropy jump in a cell.
 *
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 */
template <int dim>
double ConservationLaw<dim>::compute_max_entropy_jump(
  const Vector<double> & solution, const Cell & cell) const
{
  FEFaceValues<dim> fe_values_face(fe,
                                   face_quadrature,
                                   update_values | update_gradients |
                                     update_normal_vectors);

  std::vector<Tensor<1, dim>> gradients_face(n_q_points_face);
  std::vector<Tensor<1, dim>> gradients_face_neighbor(n_q_points_face);
  std::vector<Tensor<1, dim>> normal_vectors(n_q_points_face);
  Vector<double> entropy(n_q_points_face);

  double max_jump_in_cell = 0.0;

  // loop over faces
  for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
  {
    // determine if face is interior
    typename DoFHandler<dim>::face_iterator face = cell->face(iface);
    if (face->at_boundary() == false)
    {
      // get gradients from this cell
      fe_values_face.reinit(cell, iface);
      fe_values_face.get_function_gradients(solution, gradients_face);

      // compute entropy at each quadrature point on face
      compute_entropy(solution, fe_values_face, entropy);

      // get normal vectors
      normal_vectors = fe_values_face.get_all_normal_vectors();

      // get iterator to neighboring cell and face index of this face
      Assert(cell->neighbor(iface).state() == IteratorState::valid,
             ExcInternalError());
      typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(iface);
      const unsigned int ineighbor = cell->neighbor_of_neighbor(iface);
      Assert(ineighbor < faces_per_cell, ExcInternalError());

      // get gradients from neighboring cell
      fe_values_face.reinit(neighbor, ineighbor);
      fe_values_face.get_function_gradients(solution, gradients_face_neighbor);

      // loop over face quadrature points to determine max jump on face
      double max_jump_on_face = 0.0;
      for (unsigned int q = 0; q < n_q_points_face; ++q)
      {
        // compute difference in gradients across face
        gradients_face[q] -= gradients_face_neighbor[q];
        double jump_on_face =
          std::abs((gradients_face[q] * normal_vectors[q]) * entropy[q]);
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
void ConservationLaw<dim>::smooth_entropy_viscosity_max()
{
  // copy entropy viscosities
  CellMap entropy_viscosity_unsmoothed = entropy_viscosity_map;

  // loop over cells
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
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
        entropy_viscosity_map[cell] =
          std::max(entropy_viscosity_map[cell],
                   entropy_viscosity_unsmoothed[neighbor_cell]);
      }
    }
  }
}

/**
 * \brief Smooths the entropy visosity profile using a weighted average of the
 *        cell and its neighbors.
 */
template <int dim>
void ConservationLaw<dim>::smooth_entropy_viscosity_average()
{
  // copy entropy viscosities
  CellMap entropy_viscosity_unsmoothed = entropy_viscosity_map;

  // loop over cells
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // initialize sum for average
    double sum = parameters.entropy_viscosity_smoothing_weight *
      entropy_viscosity_unsmoothed[cell];

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
    entropy_viscosity_map[cell] =
      sum / (neighbor_count + parameters.entropy_viscosity_smoothing_weight);
  }
}

/**
 * \brief Computes low-order viscosity for each cell.
 *
 * The low-order viscosity is computed as
 * \f[
 *   \nu_K^L = c_{max} h_K \lambda_{K,max} \,,
 * \f]
 * where \f$\lambda_{K,max}\f$ is the maximum flux speed on cell \f$K\f$,
 * or if the user specifies, the low-order viscosity definition above is
 * multiplied by the local Froude number:
 * \f[
 *   \nu_K^L = c_{max} \mbox{Fr}_K h_K \lambda_{K,max} \,,
 * \f]
 * where the Froude number is taken to be the max over all quadrature points:
 * \f[
 *   \mbox{Fr}_K = \max\limits_q \frac{u_q}{a_q} \,.
 * \f]
 *
 * \param[in] using_low_order_scheme flag that the low-order viscosities are
 *            are used in a low-order scheme as opposed to an entropy
 *            viscosity scheme
 */
/*
template <int dim>
void ShallowWater<dim>::update_old_low_order_viscosity(
  const bool & using_low_order_scheme)
{
  const double c_max = this->parameters.first_order_viscosity_coef;

  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();

  if (using_low_order_scheme &&
      sw_parameters.multiply_low_order_viscosity_by_froude)
  {
    FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);

    // compute low-order viscosity for each cell
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      // extract height and momentum
      std::vector<double> height(this->n_q_points_cell);
      fe_values[height_extractor].get_function_values(this->new_solution, height);
      std::vector<Tensor<1, dim>> momentum(this->n_q_points_cell);
      fe_values[momentum_extractor].get_function_values(this->new_solution,
                                                        momentum);

      // compute fluid speed and sound speed
      std::vector<double> speed = compute_speed(height, momentum);
      std::vector<double> sound_speed = compute_sound_speed(height);

      // compute Froude number for cell by taking the maximum over
      // quadrature points
      double froude = 0.0;
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
        froude = std::max(froude, speed[q] / sound_speed[q]);

      this->first_order_viscosity_map[cell] =
        std::abs(c_max * froude * this->cell_diameter[cell] *
                 this->max_flux_speed_cell[cell]);
    }
  }
  else
  {
    // compute low-order viscosity for each cell
    for (; cell != endc; ++cell)
    {
      this->first_order_viscosity_map[cell] = std::abs(
        c_max * this->cell_diameter[cell] * this->max_flux_speed_cell[cell]);
    }
  }
}
*/

/**
 * \brief Computes entropy viscosity for each cell.
 *
 * \param[in] dt time step size
 */
template <int dim>
void ShallowWater<dim>::update_entropy_viscosities(const double & dt)
{
  // compute normalization constant for entropy viscosity
  double entropy_average = 1.0e15;
  if (sw_parameters.entropy_normalization == "average")
    entropy_average = this->compute_average_entropy(this->new_solution);

  // FE values for entropy flux
  ShallowWaterEntropyFluxFEValuesCell<dim> entropy_flux_fe_values_cell(
    this->dof_handler,
    this->triangulation,
    this->cell_quadrature,
    bathymetry_vector,
    gravity);
  ShallowWaterEntropyFluxFEValuesFace<dim> entropy_flux_fe_values_face(
    this->dof_handler,
    this->triangulation,
    this->face_quadrature,
    bathymetry_vector,
    gravity);

  // compute entropy viscosity for each cell
  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reinitialize entropy flux FE values for cell (face values will need be
    // reinitialized in compute_max_entropy_jump())
    entropy_flux_fe_values_cell.reinit(cell);

    // compute entropy normalization
    std::vector<double> entropy_normalization(this->n_q_points_cell);
    if (sw_parameters.entropy_normalization == "average")
      entropy_normalization = this->compute_entropy_normalization(
        this->new_solution, entropy_average, cell);
    else if (sw_parameters.entropy_normalization == "local")
      entropy_normalization =
        compute_local_entropy_normalization(this->new_solution, cell);
    else if (sw_parameters.entropy_normalization == "constant")
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
        entropy_normalization[q] =
          sw_parameters.constant_entropy_normalization_coefficient;
    else
      ExcNotImplemented();

    // compute entropy residual at each quadrature point on cell
    const std::vector<double> entropy_residual =
      this->compute_entropy_residual(this->new_solution,
                                     this->old_solution,
                                     entropy_flux_fe_values_cell,
                                     dt,
                                     cell);

    // compute max entropy flux jump
    const double max_entropy_jump =
      compute_max_entropy_jump(entropy_flux_fe_values_face, cell);

    // compute entropy viscosity
    double h2 = std::pow(this->cell_diameter[cell], 2);
    this->entropy_viscosity_map[cell] = 0.0;
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      this->entropy_viscosity_map[cell] = std::max(
        this->entropy_viscosity_map[cell],
        h2 * (entropy_residual_coefficient * std::abs(entropy_residual[q]) +
              jump_coefficient * max_entropy_jump) /
          entropy_normalization[q]);
  }
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
std::vector<double> EntropyViscosity<dim>::compute_entropy_residual(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const GroupFEValuesCell<dim> & entropy_flux_fe_values,
  const double & dt,
  const Cell & cell) const
{
  // FE values
  FEValues<dim> fe_values(
    *fe, *cell_quadrature, update_values | update_gradients);
  fe_values.reinit(cell);

  Vector<double> entropy_new(n_q_points_cell);
  Vector<double> entropy_old(n_q_points_cell);
  std::vector<double> entropy_residual(n_q_points_cell);

  // compute entropy of current and old solutions
  auto entropy_new = entropy->compute_entropy(new_solution, fe_values);
  auto entropy_old = entropy->compute_entropy(old_solution, fe_values);
  auto divergence_entropy_flux = entropy->compute_divergence_entropy_flux();

  // compute entropy residual at each quadrature point on cell
  for (unsigned int q = 0; q < n_q_points_cell; ++q)
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
