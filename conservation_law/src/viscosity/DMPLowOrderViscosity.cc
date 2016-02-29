/**
 * \file DMPLowOrderViscosity.cc
 * \brief Provides the function definitions for the DMPLowOrderViscosity class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
DMPLowOrderViscosity<dim>::DMPLowOrderViscosity()
{
}

/**
 * \brief Computes the maximum-principle preserving first order viscosity
 *        for each cell.
 */
template <int dim>
void ConservationLaw<dim>::update_max_principle_viscosity()
{
  // compute viscous fluxes
  compute_viscous_fluxes();

  // local dof indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells to compute first order viscosity at each quadrature point
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    first_order_viscosity_map[cell] = 0.0;
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        if (i != j)
        {
          first_order_viscosity_map[cell] = std::max(
            first_order_viscosity_map[cell],
            std::abs(viscous_fluxes(local_dof_indices[i], local_dof_indices[j])) /
              (-viscous_bilinear_forms(local_dof_indices[i],
                                       local_dof_indices[j])));
        }
      }
    }
  }
}

/**
 * \brief Computes viscous fluxes, to be used in the computation of
 *        maximum-principle preserving first order viscosity.
 *
 * Each element of the resulting matrix, \f$V_{i,j}\f$ is computed as follows:
 * \f[
 *   V_{i,j} =
 *     \int_{S_{i,j}}(\mathbf{f}'(u)\cdot\nabla\varphi_j)\varphi_i d\mathbf{x}
 * \f]
 */
template <int dim>
void ConservationLaw<dim>::compute_viscous_fluxes()
{
  viscous_fluxes = 0; // zero out matrix

  // local dof indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  FEValues<dim> fe_values(
    fe, cell_quadrature, update_values | update_gradients | update_JxW_values);
  std::vector<double> solution_values(n_q_points_cell);
  std::vector<Tensor<1, dim>> dfdu(n_q_points_cell);

  // loop over cells
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    fe_values.reinit(cell);
    fe_values.get_function_values(new_solution, solution_values);

    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // add local viscous fluxes to global matrix
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      for (int d = 0; d < dim; ++d)
        dfdu[q][d] = solution_values[q];
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          viscous_fluxes.add(local_dof_indices[i],
                             local_dof_indices[j],
                             dfdu[q] * fe_values.shape_grad(j, q) *
                               fe_values.shape_value(i, q) * fe_values.JxW(q));
        }
      }
    }
  }
}

/** \brief Gets the values and indices of nonzero elements in a sparse matrix.
 *  \param [in] matrix sparse matrix whose row will be retrieved
 *  \param [in] i index of row to be retrieved
 *  \param [out] row_values vector of values of nonzero entries of row i
 *  \param [out] row_indices vector of indices of nonzero entries of row i
 *  \param [out] n_col number of nonzero entries of row i
 */
template <int dim>
void ConservationLaw<dim>::get_matrix_row(const SparseMatrix<double> & matrix,
                                          const unsigned int & i,
                                          std::vector<double> & row_values,
                                          std::vector<unsigned int> & row_indices,
                                          unsigned int & n_col)
{
  // get first and one-past-last iterator for row
  SparseMatrix<double>::const_iterator matrix_iterator = matrix.begin(i);
  SparseMatrix<double>::const_iterator matrix_iterator_end = matrix.end(i);

  // compute number of entries in row and then allocate memory
  n_col = matrix_iterator_end - matrix_iterator;
  row_values.reserve(n_col);
  row_indices.reserve(n_col);

  // loop over columns in row
  for (; matrix_iterator != matrix_iterator_end; ++matrix_iterator)
  {
    row_values.push_back(matrix_iterator->value());
    row_indices.push_back(matrix_iterator->column());
  }
}

/**
 * \brief Checks that the local discrete max principle is satisfied.
 * \param[in] n time step index, used for reporting violations
 */
template <int dim>
bool ConservationLaw<dim>::check_DMP(const unsigned int & n) const
{
  // check that each dof value is bounded by its neighbors
  bool DMP_satisfied = true;
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // check if dof is a Dirichlet node or not
    if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) ==
        dirichlet_nodes.end())
    {
      double value_i = new_solution(i);
      if (value_i < min_values(i))
      {
        DMP_satisfied = false;
        cout << "Max principle violated at time step " << n << " with dof " << i
             << ": " << value_i << " < " << min_values(i) << std::endl;
      }
      if (value_i > max_values(i))
      {
        DMP_satisfied = false;
        cout << "Max principle violated at time step " << n << " with dof " << i
             << ": " << value_i << " > " << max_values(i) << std::endl;
      }
    }
  }

  return DMP_satisfied;
}

/** \brief Computes min and max quantities for max principle
 */
template <int dim>
void ConservationLaw<dim>::compute_max_principle_quantities()
{
  // initialize min and max values
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    max_values(i) = old_solution(i);
    min_values(i) = old_solution(i);
  }

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // find min and max values on cell - start by initializing to arbitrary cell
    double max_cell = old_solution(local_dof_indices[0]);
    double min_cell = old_solution(local_dof_indices[0]);
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      double value_j = old_solution(local_dof_indices[j]);
      max_cell = std::max(max_cell, value_j);
      min_cell = std::min(min_cell, value_j);
    }

    // update the max and min values of neighborhood of each dof
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      unsigned int i = local_dof_indices[j]; // global index
      max_values(i) = std::max(max_values(i), max_cell);
      min_values(i) = std::min(min_values(i), min_cell);
    }
  }
}

/**
 * \brief Gets a list of dofs subject to Dirichlet boundary conditions.
 */
template <int dim>
void ConservationLaw<dim>::get_dirichlet_nodes()
{
  // get map of Dirichlet dof indices to Dirichlet values
  std::map<unsigned int, double> boundary_values;

  // clear Dirichlet nodes vector
  dirichlet_nodes.clear();

  if (boundary_conditions_type == "dirichlet")
    // loop over boundary IDs
    for (unsigned int boundary = 0; boundary < n_boundaries; ++boundary)
      // loop over components
      for (unsigned int component = 0; component < n_components; ++component)
      {
        // create mask to prevent function from being applied to all components
        std::vector<bool> component_mask(n_components, false);
        component_mask[component] = true;

        // get boundary values map using interpolate_boundary_values function
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 boundary,
                                                 ZeroFunction<dim>(n_components),
                                                 boundary_values,
                                                 component_mask);

        // extract dof indices from map
        for (std::map<unsigned int, double>::iterator it =
               boundary_values.begin();
             it != boundary_values.end();
             ++it)
          dirichlet_nodes.push_back(it->first);
      }
}
