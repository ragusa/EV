/**
 * \file FunctionDoFBounds.cc
 * \brief Provides the function definitions for the FunctionDoFBounds class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] function_  function
 * \param[in] function_is_time_dependent_  flag that function is time-dependent
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] cell_quadrature_  cell quadrature
 */
template <int dim>
FunctionDoFBounds<dim>::FunctionDoFBounds(
  Function<dim> & function_,
  const bool & function_is_time_dependent_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const QGauss<dim> & cell_quadrature_)
  : DoFBounds<dim>(dof_handler_, fe_),
    function(&function_),
    function_is_time_dependent(function_is_time_dependent_),
    cell_quadrature(&cell_quadrature_),
    n_q_points_cell(cell_quadrature_.size())
{
}

/**
 * \brief Sets time for function if applicable and computes function bounds.
 *
 * \param[in] time  time at which to evaluate function if it is time-dependent
 */
template <int dim>
void FunctionDoFBounds<dim>::update(const double & time)
{
  // set time of function if applicable
  if (function_is_time_dependent)
    function->set_time(time);

  // update function bounds
  compute_function_bounds();
}

/**
 * \brief Computes the minimum and maximum of function
 *        in the neighborhood of each degree of freedom.
 */
template <int dim>
void FunctionDoFBounds<dim>::compute_function_bounds()
{
  // initialize min and max values
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    this->lower[i] = 1.0e15;
    this->upper[i] = -1.0e15;
  }

  // FE values
  FEValues<dim> fe_values(*this->fe, *cell_quadrature, update_quadrature_points);

  // loop over cells
  Cell cell = this->dof_handler->begin_active(), endc = this->dof_handler->end();
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);
  for (; cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // get quadrature points on cell
    std::vector<Point<dim>> points(n_q_points_cell);
    points = fe_values.get_quadrature_points();

    // find min and max values on cell
    double min_cell = 1.0e15;
    double max_cell = -1.0e15;
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      // update min and max values on cell
      const double value_q = function->value(points[q]);
      min_cell = std::min(min_cell, value_q);
      max_cell = std::max(max_cell, value_q);
    }

    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // update min and max values in neighborhood of each dof
    for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
    {
      unsigned int i = local_dof_indices[j]; // global index
      this->lower[i] = std::min(this->lower[i], min_cell);
      this->upper[i] = std::max(this->upper[i], max_cell);
    }
  }
}
