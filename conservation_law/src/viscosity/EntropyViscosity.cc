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
 * \param[in] use_in_laplacian_term_  flag that viscosity is to be used
 *            in a Laplacian diffusion term
 */
template <int dim>
EntropyViscosity<dim>::EntropyViscosity(
  const RunParameters<dim> & parameters_,
  const std::shared_ptr<Entropy<dim>> & entropy_,
  const CellMap & cell_diameter_,
  const FESystem<dim> & fe_,
  const DoFHandler<dim> & dof_handler_,
  const QGauss<dim> & cell_quadrature_,
  const QGauss<dim - 1> & face_quadrature_,
  const bool & use_in_laplacian_term_)
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
  // point viscosity multiplier function pointer to appropriate function
  if (use_in_laplacian_term_)
    compute_viscosity_multiplier_ptr =
      &EntropyViscosity<dim>::compute_viscosity_multiplier_laplacian;
  else
    compute_viscosity_multiplier_ptr =
      &EntropyViscosity<dim>::compute_viscosity_multiplier_graphtheoretic;
}

/**
 * \brief Computes viscosity values.
 *
 * When using Laplacian diffusion, the entropy viscosity for cell \f$K\f$
 * at time \f$n\f$ is computed as
 * \f[
 *   \nu^{\eta,n}_K \equiv h_K^2 \max\limits_{q\in Q(K)}\left(
 *     \frac{c_R |R_q^n| + c_J J_K^n}{\hat{\eta}_q}
 *     \right)\,.
 * \f]
 * When using graph-theoretic diffusion, the leading \f$h_K^2\f$ coefficient
 * is dropped:
 * \f[
 *   \nu^{\eta,n}_K \equiv \max\limits_{q\in Q(K)}\left(
 *     \frac{c_R |R_q^n| + c_J J_K^n}{\hat{\eta}_q}
 *     \right)\,.
 * \f]
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
    // const auto entropy_residual =
    //  compute_entropy_residual(new_solution, old_solution, dt, cell);
    const std::vector<double> entropy_residual =
      entropy->compute_entropy_residual(new_solution, old_solution, dt, cell);

    // compute maximum entropy residual
    double max_entropy_residual = 0.0;
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      max_entropy_residual =
        std::max(max_entropy_residual, std::max(0.0, entropy_residual[q]));

    // compute max entropy gradient jump in cell
    // const double max_entropy_jump = compute_max_entropy_jump(cell);
    const double max_entropy_jump =
      entropy->compute_max_entropy_jump(new_solution, cell);

    // compute normalization at each quadrature point in cell
    const std::vector<double> entropy_normalization =
      entropy->compute_entropy_normalization(new_solution, cell);

    // compute entropy viscosity
    double viscosity_cell = 0.0;
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      viscosity_cell = std::max(viscosity_cell,
                                (residual_coefficient * max_entropy_residual +
                                 jump_coefficient * max_entropy_jump) /
                                  entropy_normalization[q]);
    }

    // compute multiplier
    const double multiplier = (this->*compute_viscosity_multiplier_ptr)(cell);

    // store value
    this->values[cell] = multiplier * viscosity_cell;
  }
}

/**
 * \brief Computes viscosity multiplier for use in Laplacian diffusion term.
 *
 * \param[in] cell  cell iterator
 *
 * \return viscosity multiplier for cell
 */
template <int dim>
double EntropyViscosity<dim>::compute_viscosity_multiplier_laplacian(
  const Cell & cell) const
{
  return std::pow(cell_diameter->at(cell), 2);
}

/**
 * \brief Computes viscosity multiplier for use in graph-theoretic diffusion term.
 *
 * \param[in] cell  cell iterator
 *
 * \return viscosity multiplier for cell
 */
template <int dim>
double EntropyViscosity<dim>::compute_viscosity_multiplier_graphtheoretic(
  const Cell &) const
{
  return 1.0;
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
