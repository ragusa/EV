/**
 * \file GroupFEValuesCell.cc
 * \brief Provides the function definitions for the GroupFEValuesCell class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_components_solution_ number of solution components
 * \param[in] n_components_function_ number of function components
 * \param[in] solution_dof_handler_ DoF handler for solution
 * \param[in] triangulation_ triangulation
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] solution_ solution vector
 * \param[in] aux_vector_ optional auxiliary vector
 */
template <int dim, typename FunctionType>
GroupFEValuesCell<dim,FunctionType>::GroupFEValuesCell(
  const unsigned int & n_components_solution_,
  const unsigned int & n_components_function_,
  const DoFHandler<dim> & solution_dof_handler_,
  const Triangulation<dim> & triangulation_,
  const QGauss<dim> & cell_quadrature_,
  const Vector<double> & solution_,
  const Vector<double> & aux_vector_)
  : GroupFEValuesBase<dim,FunctionType>(n_components_solution_,
                           n_components_function_,
                           solution_dof_handler_,
                           triangulation_,
                           solution_,
                           aux_vector_),
    cell_quadrature(cell_quadrature_),
    n_quadrature_points(cell_quadrature.size()),
    function_fe_values_cell(
      this->fe, cell_quadrature, update_values | update_gradients)
{
}

/**
 * \brief Reinitializes function FE values for cell.
 *
 * \param[in] solution_cell solution cell iterator
 */
template <int dim, typename FunctionType>
void GroupFEValuesCell<dim,FunctionType>::reinit(const cell_iterator & solution_cell)
{
  // reinitialize function FE values with the function cell corresponding
  // to the solution cell
  function_fe_values_cell.reinit(
    this->solution_cell_to_function_cell_map[solution_cell]);
}

/**
 * \brief Computes function values at quadrature points in a cell.
 *
 * \param[out] function_values vector of function values at quadrature points
 *             in a cell
 */
template <int dim, typename FunctionType>
void GroupFEValuesCell<dim,FunctionType>::get_function_values(
  std::vector<double> & function_values) const
{
  function_fe_values_cell[this->function_extractor].get_function_values(
    this->function_dof_values, function_values);
}
