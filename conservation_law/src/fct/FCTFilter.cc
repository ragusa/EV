/**
 * \file FCTFilter.cc
 * \brief Provides the function definitions for the FCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 */
template <int dim>
FCTFilter<dim>::FCTFilter(
  const RunParameters & run_parameters_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const std::map<unsigned int, double> & dirichlet_values_)
  : solution_bounds(dof_handler_, fe_),
    antidiffusion_bounds(dof_handler_, fe_),
    limiter(limiter_),
    dof_handler(&dof_handler_),
    n_dofs(dof_handler_.n_dofs()),
    n_components(dof_handler_.get_fe().n_components()),
    dofs_per_cell(dof_handler_.get_fe().dofs_per_cell),
    dofs_per_cell_per_component(dofs_per_cell / n_components),
    do_enforce_antidiffusion_bounds_signs(
      run_parameters_.enforce_antidiffusion_bounds_signs),
    dirichlet_values(&dirichlet_values_)
{
}

/**
 * \brief Checks to see if the FCT bounds were satisfied.
 *
 * \param[in] new_solution  new solution vector \f$\mathbf{U}^{n+1}\f$
 *
 * \return flag that FCT bounds were satisfied for all filters
 */
template <int dim>
bool FCTFilter<dim>::check_bounds(const Vector<double> & new_solution)
{
  return solution_bounds.check_bounds(new_solution);
}

/**
 * \brief Gets the lower solution bound.
 *
 * \return lower solution bound vector
 */
template <int dim>
Vector<double> FCTFilter<dim>::get_lower_solution_bound() const
{
  return solution_bounds.lower;
}

/**
 * \brief Gets the upper solution bound.
 *
 * \return upper solution bound vector
 */
template <int dim>
Vector<double> FCTFilter<dim>::get_upper_solution_bound() const
{
  return solution_bounds.upper;
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

/**
 * \brief Enforces the necessary signs of the antidiffusion bounds.
 *
 * The antidiffusion bounds should have the following properties for use
 * in the FCT algorithm:
 * \f[
 *   Q_i^- \leq 0 \,,
 * \f]
 * \f[
 *   Q_i^+ \geq 0 \,.
 * \f]
 * This function enforces these properties by performing the following
 * operations:
 * \f[
 *   Q_i^- \gets \min(0, Q_i^-) \,,
 * \f]
 * \f[
 *   Q_i^+ \gets \max(0, Q_i^+) \,.
 * \f]
 */
template <int dim>
void FCTFilter<dim>::enforce_antidiffusion_bounds_signs()
{
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    antidiffusion_bounds.lower[i] = std::min(0.0, antidiffusion_bounds.lower[i]);
    antidiffusion_bounds.upper[i] = std::max(0.0, antidiffusion_bounds.upper[i]);
  }
}

/**
 * \briefs Checks that the signs of the antidiffusion bounds \f$Q_i^\pm\f$ are
 *         correct.
 */
template <int dim>
void FCTFilter<dim>::check_antidiffusion_bounds_signs() const
{
#ifdef DEBUG
  // small number for checking sign
  const double small = 1.0e-15;
#endif

  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // if not a Dirichlet node
    if (dirichlet_values->find(i) == dirichlet_values->end())
    {
      Assert(
        antidiffusion_bounds.lower[i] < small,
        ExcAntidiffusionLowerBoundPositive(i, antidiffusion_bounds.lower[i]));
      Assert(
        antidiffusion_bounds.upper[i] > -small,
        ExcAntidiffusionUpperBoundNegative(i, antidiffusion_bounds.upper[i]));
    }
  }
}
