/**
 * \file FCTFilter.cc
 * \brief Provides the function definitions for the FCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 */
template <int dim>
FCTFilter<dim>::FCTFilter(const std::shared_ptr<Limiter> limiter_,
                          const DoFHandler<dim> & dof_handler_)
  : solution_bounds(dof_handler_.n_dofs()),
    antidiffusion_bounds(dof_handler_.n_dofs()),
    limiter(limiter_),
    dof_handler(&dof_handler_),
    n_dofs(dof_handler_.n_dofs()),
    n_components(dof_handler_.get_fe().n_components()),
    dofs_per_cell(dof_handler_.get_fe().dofs_per_cell),
    dofs_per_cell_per_component(dofs_per_cell / n_components)
{
}

/**
 * \brief Computes the minimum and maximum degree of freedom values in the
 *        support of each degree of freedom.
 *
 * \param[in] dof_vector  degree of freedom vector \f$\mathbf{U}\f$
 * \param[out] min_values  minimum values of solution:
 *             \f$U_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} U_j\f$
 * \param[out] max_values  maximum values of solution:
 *             \f$U_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} U_j\f$
 */
template <int dim>
void FCTFilter<dim>::compute_min_and_max_of_dof_vector(
  const Vector<double> & dof_vector,
  Vector<double> & min_values,
  Vector<double> & max_values) const
{
  // initialize solution min and max
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    min_values(i) = dof_vector(i);
    max_values(i) = dof_vector(i);
  }

  // loop over cells
  Cell cell = dof_handler->begin_active(), endc = dof_handler->end();
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // loop over dof_vector components
    for (unsigned int m = 0; m < n_components; ++m)
    {
      // find min and max values on cell for this component
      double max_cell = dof_vector(local_dof_indices[m]);
      double min_cell = dof_vector(local_dof_indices[m]);
      for (unsigned int j = 0; j < dofs_per_cell_per_component; ++j)
      {
        // consider old states
        const double value_j =
          dof_vector(local_dof_indices[j * n_components + m]);
        max_cell = std::max(max_cell, value_j);
        min_cell = std::min(min_cell, value_j);
      }

      // update the max and min values of neighborhood of each dof
      for (unsigned int j = 0; j < dofs_per_cell_per_component; ++j)
      {
        unsigned int i = local_dof_indices[j * n_components + m]; // global index
        max_values(i) = std::max(max_values(i), max_cell);
        min_values(i) = std::min(min_values(i), min_cell);
      }
    }
  }
}
