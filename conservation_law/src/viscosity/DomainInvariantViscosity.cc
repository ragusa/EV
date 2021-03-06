/**
 * \file DomainInvariantViscosity.cc
 * \brief Provides the function definitions for the DomainInvariantViscosity
 *        class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] max_wave_speed_ max wave speed object
 * \param[in] gradient_matrix_ gradient matrix object
 * \param[in] fe_ finite element system
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] triangulation_ triangulation
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] n_components_ number of solution components
 * \param[in] viscosity_multiplier_ viscosity multiplier
 */
template <int dim>
DomainInvariantViscosity<dim>::DomainInvariantViscosity(
  const std::shared_ptr<MaxWaveSpeed<dim>> & max_wave_speed_,
  const std::shared_ptr<GradientMatrix<dim>> & gradient_matrix_,
  const FESystem<dim> & fe_,
  const DoFHandler<dim> & dof_handler_,
  const Triangulation<dim> & triangulation_,
  const QGauss<dim> & cell_quadrature_,
  const unsigned int & n_components_,
  const std::shared_ptr<ViscosityMultiplier<dim>> & viscosity_multiplier_)
  : Viscosity<dim>(dof_handler_),
    max_wave_speed(max_wave_speed_),
    gradient_matrix(gradient_matrix_),
    fe_scalar(1),
    fe(&fe_),
    dof_handler_scalar(triangulation_),
    dof_handler(&dof_handler_),
    cell_quadrature(&cell_quadrature_),
    n_components(n_components_),
    n_dofs_scalar(dof_handler->n_dofs() / n_components),
    dofs_per_cell_scalar(fe_scalar.dofs_per_cell),
    viscosity_multiplier(viscosity_multiplier_)
{
  // reinitialize
  reinitialize();
}

/**
 * \brief Reinitializes.
 */
template <int dim>
void DomainInvariantViscosity<dim>::reinitialize()
{
  // distribute degrees of freedom for scalar case
  dof_handler_scalar.clear();
  dof_handler_scalar.distribute_dofs(fe_scalar);

  // create sparsity pattern and reinitialize sparse matrices
  DynamicSparsityPattern dsp(n_dofs_scalar);
  DoFTools::make_sparsity_pattern(dof_handler_scalar, dsp);
  sparsity.copy_from(dsp);
  viscous_sums.reinit(sparsity);

  // compute viscous bilinear form sum matrix
  compute_graph_theoretic_sums();
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
void DomainInvariantViscosity<dim>::update(const Vector<double> & new_solution,
                                           const Vector<double> &,
                                           const double &,
                                           const unsigned int &)
{
  // FE values
  FEValues<dim> fe_values(*fe, *cell_quadrature, update_values);

  // global degree of freedom indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell_scalar);

  // loop over cells
  Cell cell_scalar = dof_handler_scalar.begin_active();
  Cell endc_scalar = dof_handler_scalar.end();
  Cell cell = this->dof_handler->begin_active();
  for (; cell_scalar != endc_scalar; ++cell_scalar, ++cell)
  {
    // reinitialize fe values for cell
    fe_values.reinit(cell);

    // get multiplier
    const double multiplier =
      viscosity_multiplier->get_multiplier(fe_values, new_solution);

    // get global degree of freedom indices
    cell_scalar->get_dof_indices(local_dof_indices);

    // reset viscosity for max() function
    this->values[cell] = 0.0;

    // loop over pairs of local degrees of freedom
    for (unsigned int i_local = 0; i_local < dofs_per_cell_scalar; ++i_local)
      for (unsigned int j_local = 0; j_local < dofs_per_cell_scalar; ++j_local)
        if (j_local != i_local)
        {
          // get global scalar indices for i and j
          const unsigned int i = local_dof_indices[i_local];
          const unsigned int j = local_dof_indices[j_local];

          // extract solution
          std::vector<double> solution_i(n_components);
          std::vector<double> solution_j(n_components);
          for (unsigned int m = 0; m < n_components; ++m)
          {
            solution_i[m] = new_solution[i * n_components + m];
            solution_j[m] = new_solution[j * n_components + m];
          }

          // get normal vector
          auto normal = gradient_matrix->get_normal(i, j);

          // compute maximum wave speed
          const double max_wave_speed_value =
            max_wave_speed->compute(solution_i, solution_j, normal);

          // gradient matrix entry
          const double gradient_norm = gradient_matrix->get_gradient_norm(i, j);

          // compute viscosity
          const double viscosity = -(max_wave_speed_value * gradient_norm) /
            viscous_sums(i, j) * multiplier;

          // update max viscosity for this cell
          this->values[cell] = std::max(this->values[cell], viscosity);
        }
  }
}

/**
 * \brief Computes sums of the graph-theoretic viscous bilinear forms.
 *
 * Each element of the resulting matrix, \f$B_{i,j}\f$ is computed as follows:
 * \f[
 *   B_{i,j} = \sum_{K:D_K\subset S_{i,j}}b_K(\varphi_i,\varphi_j)
 * \f]
 */
template <int dim>
void DomainInvariantViscosity<dim>::compute_graph_theoretic_sums()
{
  // reset matrix to zero
  viscous_sums = 0;

  // local dof indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell_scalar);

  // loop over cells
  Cell cell = dof_handler_scalar.begin_active(), endc = dof_handler_scalar.end();
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // get cell volume
    const double cell_volume = cell->measure();

    // add local bilinear forms to global matrix
    for (unsigned int i = 0; i < dofs_per_cell_scalar; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell_scalar; ++j)
      {
        // compute local viscous bilinear form
        double b_cell = 0.0;
        if (j == i)
          b_cell = cell_volume;
        else
          b_cell = -1.0 / (dofs_per_cell_scalar - 1.0) * cell_volume;

        // add local contribution to sum
        viscous_sums.add(local_dof_indices[i], local_dof_indices[j], b_cell);
      }
    }
  }
}
