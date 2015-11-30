/**
 * \file DomainInvariantViscosity.cc
 * \brief Provides the function definitions for the DomainInvariantViscosity
 *        class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
DomainInvariantViscosity<dim>::DomainInvariantViscosity(
  const Triangulation<dim> & triangulation)
  : fe(1),
    dof_handler(triangulation_),
    triangulation(&triangulation_),
    lines_per_cell(GeometryInfo<dim>::lines_per_cell)
{
  // reinitialize here?
}

/**
 * \brief Reinitializes.
 */
template <int dim>
void DomainInvariantViscosity<dim>::reinitialize()
{
  // distribute DoFs
  dof_handler.clear();
  dof_handler.distribute_dofs(fe);

  // create sparsity pattern and reinitialize sparse matrices
  DynamicSparsityPattern dsp(n_dofs);
  DoFTools::make_sparsity_pattern(dof_handler.n_dofs(), dsp);
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
                                           const Vector<double> & old_solution,
                                           const double & dt,
                                           const unsigned int &)
{
  // currently this is only valid for 1-D because in 2-D a line may
  // correspond to more than 1 cell
  Assert(dim == 1, ExcInvalidInDimension(dim))

    // dof indices for a line
    std::vector<unsigned int> local_dof_indices(2);

  // clear max viscosity map
  max_viscosity.clear();

  // loop over cells
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // loop over lines of cell
    for (unsigned int line = 0; line < lines_per_cell; ++line)
    {
      // get line iterator
      const Line line_iterator = cell->line(line);

      // determine if line iterator exists in map; if not, need to
      // compute viscosity
      if (max_viscosity.find(line_iterator) == max_viscosity.end())
      {
        // get dof indices on line
        line_iterator->get_dof_indices(local_dof_indices);
        const unsigned int i = local_dof_indices[0];
        const unsigned int j = local_dof_indices[1];

        // extract solution
        std::vector<double> solution_i(n_components);
        std::vector<double> solution_j(n_components);
        for (unsigned int m = 0; m < n_components; ++m)
        {
          solution_i[m] = new_solution[i * n_components + m];
          solution_j[m] = new_solution[j * n_components + m];
        }

        // compute normal vectors
        const Tensor<1, dim> normal_ij[i_line] =
          c_ij[i_line] / c_ij[i_line].norm();
        const Tensor<1, dim> normal_ji[i_line] =
          c_ji[i_line] / c_ji[i_line].norm();

        // compute maximum wave speeds
        const double max_wave_speed_ij =
          max_wave_speed->compute(normal_ij[i_line], solution_i, solution_j);
        const double max_wave_speed_ji =
          max_wave_speed->compute(normal_ji[i_line], solution_j, solution_i);

        // compute viscosities
        const double viscosity_ij =
          -(max_wave_speed_ij * c_ij_norm[i_line]) / viscous_sums(i, j);
        const double viscosity_ji =
          -(max_wave_speed_ji * c_ji_norm[i_line]) / viscous_sums(j, i);

        // store max viscosity in line map
        max_viscosity[line_iterator] = std::max(viscosity_ij, viscosity_ji);
      }

      // update max viscosity for this cell
      this->values[cell] =
        std::max(this->values[cell], max_viscosity[line_iterator]);
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
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // get cell volume
    const double cell_volume = cell->measure();

    // add local bilinear forms to global matrix
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        // compute local viscous bilinear form
        double b_cell = 0.0;
        if (j == i)
          b_cell = cell_volume;
        else
          b_cell = -1.0 / (dofs_per_cell - 1.0) * cell_volume;

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
  FEValues<dim> fe_values(
    fe, cell_quadrature, update_values | update_gradients | update_JxW_values);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reinitialize fe values for cell
    fe_values.reinit(cell);

    // get DoF indices
    cell->get_dof_indices(local_dof_indices);

    // loop over dimension
    for (unsigned int d = 0; d < dim; ++d)
      // loop over test function i
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        // loop over test function j
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
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

  // compute normals of each gradient entry
  SparsityPattern::const_iterator it = sparsity.begin(), it_end = sparsity.end();
  for (; it != it_end; ++it)
  {
  }
}
