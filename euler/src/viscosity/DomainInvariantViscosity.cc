/**
 * \file DomainInvariantViscosity.cc
 * \brief Provides the function definitions for the DomainInvariantViscosity
 *        class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] max_wave_speed_ max wave speed object
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] triangulation_ triangulation
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] n_components_ number of solution components
 */
template <int dim>
DomainInvariantViscosity<dim>::DomainInvariantViscosity(
  const std::shared_ptr<MaxWaveSpeed<dim>> & max_wave_speed_,
  const DoFHandler<dim> & dof_handler_,
  const Triangulation<dim> & triangulation_,
  const QGauss<dim> & cell_quadrature_,
  const unsigned int & n_components_)
  : Viscosity<dim>(dof_handler_),
    max_wave_speed(max_wave_speed_),
    fe(1),
    dof_handler_scalar(triangulation_),
    dof_handler(&dof_handler_),
    cell_quadrature(&cell_quadrature_),
    n_components(n_components_),
    n_dofs_scalar(dof_handler->n_dofs() / n_components),
    dofs_per_cell_scalar(fe.dofs_per_cell),
    n_q_points_cell(cell_quadrature_.size())
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
  dof_handler_scalar.distribute_dofs(fe);

  // create sparsity pattern and reinitialize sparse matrices
  DynamicSparsityPattern dsp(n_dofs_scalar);
  DoFTools::make_sparsity_pattern(dof_handler_scalar, dsp);
  sparsity.copy_from(dsp);
  viscous_sums.reinit(sparsity);
  for (unsigned int d = 0; d < dim; ++d)
    gradients[d].reinit(sparsity);
  gradient_norms.reinit(sparsity);

  // compute viscous bilinear form sum matrix
  compute_graph_theoretic_sums();

  // compute gradients
  compute_gradients_and_normals();
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
  // global degree of freedom indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell_scalar);

  // loop over cells
  Cell cell_scalar = dof_handler_scalar.begin_active();
  Cell endc_scalar = dof_handler_scalar.end();
  Cell cell = this->dof_handler->begin_active();
  for (; cell_scalar != endc_scalar; ++cell_scalar, ++cell)
  {
    // get global degree of freedom indices
    cell_scalar->get_dof_indices(local_dof_indices);

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

          // create tensor for normal vector
          Tensor<1, dim> normal;
          for (unsigned int d = 0; d < dim; ++d)
            normal[d] = gradients[d](i, j);

          // compute maximum wave speed
          const double max_wave_speed_value =
            max_wave_speed->compute(solution_i, solution_j, normal);

          // compute viscosity
          const double viscosity =
            -(max_wave_speed_value * gradient_norms(i, j)) / viscous_sums(i, j);

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

/**
 * \brief Computes gradients and normal vectors along each edge.
 */
template <int dim>
void DomainInvariantViscosity<dim>::compute_gradients_and_normals()
{
  // reset matrices
  for (unsigned int d = 0; d < dim; ++d)
    gradients[d] = 0;
  gradient_norms = 0;

  FEValues<dim> fe_values(
    fe, *cell_quadrature, update_values | update_gradients | update_JxW_values);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell_scalar);

  // loop over cells
  Cell cell = dof_handler_scalar.begin_active(), endc = dof_handler_scalar.end();
  for (; cell != endc; ++cell)
  {
    // reinitialize fe values for cell
    fe_values.reinit(cell);

    // get DoF indices
    cell->get_dof_indices(local_dof_indices);

    // loop over dimension
    for (unsigned int d = 0; d < dim; ++d)
      // loop over test function i
      for (unsigned int i = 0; i < dofs_per_cell_scalar; ++i)
        // loop over test function j
        for (unsigned int j = 0; j < dofs_per_cell_scalar; ++j)
        {
          // initialize integral to zero
          double entry_cell_i_j = 0.0;

          // loop over quadrature points
          for (unsigned int q = 0; q < n_q_points_cell; ++q)
            entry_cell_i_j += fe_values.shape_value(i, q) *
              fe_values.shape_grad(j, q)[d] * fe_values.JxW(q);

          // add contribution to global matrix
          gradients[d].add(
            local_dof_indices[i], local_dof_indices[j], entry_cell_i_j);
        }
  }

  // initialize matrix iterators
  SparseMatrix<double>::iterator it = gradient_norms.begin();
  SparseMatrix<double>::iterator it_end = gradient_norms.end();
  std::vector<SparseMatrix<double>::iterator> it_grad;
  for (unsigned int d = 0; d < dim; ++d)
    it_grad.push_back(gradients[d].begin());

  for (; it != it_end; ++it)
  {
    // compute gradient norms
    for (unsigned int d = 0; d < dim; ++d)
      it->value() += it_grad[d]->value() * it_grad[d]->value();
    it->value() = std::sqrt(it->value());

    // compute normalized gradient vector
    for (unsigned int d = 0; d < dim; ++d)
    {
      // normalize this component of gradient
      it_grad[d]->value() = it_grad[d]->value() / it->value();

      // increment iterator
      ++it_grad[d];
    }
  }
}
