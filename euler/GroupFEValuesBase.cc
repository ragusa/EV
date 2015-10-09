/**
 * \file GroupFEValuesBase.cc
 * \brief Provides the function definitions for the GroupFEValuesBase class.
 */

/**
 * \brief Constructor.
 *
 * In addition to construction, this function distributes DoFs and maps
 * solution cell iterators to function cell iterators. Derived classes must
 * compute the function DoF values (it is convenient to do this in their
 * constructors) by calling compute_function_dof_values(). This cannot
 * be called in the base class constructor because it would be called
 * before the function is defined, which occurs later in the construction
 * of the derived class.
 *
 * \param[in] n_components_solution_ number of solution components
 * \param[in] n_components_function_ number of function components
 * \param[in] solution_dof_handler_ DoF handler for solution
 * \param[in] triangulation_ triangulation
 * \param[in] solution_ solution vector
 * \param[in] aux_vector_ optional auxiliary vector
 */
template <int dim, typename FunctionType>
GroupFEValuesBase<dim,FunctionType>::GroupFEValuesBase(
  const unsigned int & n_components_solution_,
  const unsigned int & n_components_function_,
  const DoFHandler<dim> & solution_dof_handler_,
  const Triangulation<dim> & triangulation_,
  const Vector<double> & solution_,
  const Vector<double> & aux_vector_) :
  n_components_solution(n_components_solution_),
  n_components_function(n_components_function_),
  fe(FE_Q<dim>(1), n_components_function),
  function_extractor(0),
  dof_handler(triangulation_),
  n_function_dofs_per_cell(fe.dofs_per_cell),
  n_solution_dofs_per_cell(n_function_dofs_per_cell*n_components_solution),
  solution(&solution_),
  aux_vector(&aux_vector_)
{
  // distribute DoFs and initialize cell iterators
  dof_handler.distribute_dofs(fe);
  n_dofs = dof_handler.n_dofs();

  // map solution cell iterator to function cell iterator
  cell_iterator solution_cell = solution_dof_handler_.begin_active(),
    solution_endc = solution_dof_handler_.end(),
    function_cell = dof_handler.begin_active();
  for (; solution_cell != solution_endc; ++solution_cell, ++function_cell)
    solution_cell_to_function_cell_map[solution_cell] = function_cell;

  // resize function DoFs vector
  function_dof_values.reinit(n_dofs);
}

/**
 * \brief Computes the function DoF values.
 *
 * This function must be called by the derived class (it is best to do so in
 * the constructor).
 */
template <int dim, typename FunctionType>
void GroupFEValuesBase<dim,FunctionType>::compute_function_dof_values()
{
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // extract all components of solution at DoF support point
    std::vector<double> solution_at_support_point(n_components_solution);
    for (unsigned int j = 0; j < n_components_solution; ++j)
      solution_at_support_point[j] = (*solution)[i*n_components_solution + j];

    // compute function at DoF support point
    std::vector<double> function_at_support_point = function(solution_at_support_point);
    for (unsigned int j = 0; j < n_components_function; ++j)
      function_dof_values[i*n_components_function + j] = function_at_support_point[j];
  }
}

/**
 * \brief Gets the function DoF values.
 *
 * \return function values at DoF support points
 */
template <int dim, typename FunctionType>
Vector<double> GroupFEValuesBase<dim,FunctionType>::get_function_dof_values() const
{
  return function_dof_values;
}
