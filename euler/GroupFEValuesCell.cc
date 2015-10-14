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
template <int dim, bool is_scalar>
GroupFEValuesCell<dim, is_scalar>::GroupFEValuesCell(
  const unsigned int & n_components_solution_,
  const unsigned int & n_components_function_,
  const DoFHandler<dim> & solution_dof_handler_,
  const Triangulation<dim> & triangulation_,
  const QGauss<dim> & cell_quadrature_,
  const Vector<double> & solution_,
  const Vector<double> & aux_vector_)
  : GroupFEValuesBase<dim, is_scalar>(n_components_solution_,
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
template <int dim, bool is_scalar>
void GroupFEValuesCell<dim, is_scalar>::reinit(const Cell & solution_cell)
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
template <int dim, bool is_scalar>
void GroupFEValuesCell<dim, is_scalar>::get_function_values(
  std::vector<ValueType> & function_values) const
{
  function_fe_values_cell[this->function_extractor].get_function_values(
    this->function_dof_values, function_values);
}

/**
 * \brief Computes function gradients at quadrature points in a cell.
 *
 * \param[out] function_gradients vector of function gradients at quadrature
 *             points in a cell
 */
template <int dim, bool is_scalar>
void GroupFEValuesCell<dim, is_scalar>::get_function_gradients(
  std::vector<GradientType> & function_gradients) const
{
  function_fe_values_cell[this->function_extractor].get_function_gradients(
    this->function_dof_values, function_gradients);
}

/**
 * \brief Computes function divergences at quadrature points in a cell.
 *
 * \return vector of function divergences at quadrature points in a cell
 */
template <int dim, bool is_scalar>
std::vector<double> GroupFEValuesCell<dim, is_scalar>::get_function_divergences() const
{
  std::vector<double> function_divergences(n_quadrature_points);
  function_fe_values_cell[this->function_extractor].get_function_divergences(
    this->function_dof_values, function_divergences);

  return function_divergences;
}
