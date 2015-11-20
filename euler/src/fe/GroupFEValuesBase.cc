/**
 * \file GroupFEValuesBase.cc
 * \brief Provides the function definitions for the GroupFEValuesBase class.
 */

/**
 * \brief Constructor.
 *
 * In addition to construction, this function distributes DoFs and maps
 * solution cell iterators to function cell iterators.
 *
 * \param[in] n_components_solution_ number of solution components
 * \param[in] n_components_function_ number of function components
 * \param[in] solution_dof_handler_ DoF handler for solution
 * \param[in] triangulation_ triangulation
 * \param[in] aux_vector_ optional auxiliary vector
 */
template <int dim, bool is_scalar>
GroupFEValuesBase<dim, is_scalar>::GroupFEValuesBase(
  const unsigned int & n_components_solution_,
  const unsigned int & n_components_function_,
  const DoFHandler<dim> & solution_dof_handler_,
  const Triangulation<dim> & triangulation_,
  const Vector<double> & aux_vector_)
  : n_components_solution(n_components_solution_),
    n_components_function(n_components_function_),
    fe(FE_Q<dim>(1), n_components_function),
    function_extractor(0),
    dof_handler(triangulation_),
    n_function_dofs_per_cell(fe.dofs_per_cell),
    n_solution_dofs_per_cell(n_function_dofs_per_cell * n_components_solution),
    aux_vector(&aux_vector_),
    function_values_initialized(false)
{
  // distribute DoFs and initialize cell iterators
  dof_handler.distribute_dofs(fe);
  n_dofs = dof_handler.n_dofs();

  // compute number of support points
  n_support_points = n_dofs / n_components_function;

  // map solution cell iterator to function cell iterator
  Cell solution_cell = solution_dof_handler_.begin_active(),
       solution_endc = solution_dof_handler_.end(),
       function_cell = dof_handler.begin_active();
  for (; solution_cell != solution_endc; ++solution_cell, ++function_cell)
    solution_cell_to_function_cell_map[solution_cell] = function_cell;

  // resize function DoFs vector
  function_dof_values.reinit(n_dofs);
}

/**
 * \brief Computes the function DoF values. This function must be called
 *        before any get-values functions are called.
 *
 * \param[in] solution solution vector
 */
template <int dim, bool is_scalar>
void GroupFEValuesBase<dim, is_scalar>::reinitialize(
  const Vector<double> & solution)
{
  // set flag to signal that function DoF values have been initialized
  function_values_initialized = true;

  // loop over support points
  for (unsigned int i = 0; i < n_support_points; ++i)
  {
    // extract all components of solution at DoF support point
    std::vector<double> solution_at_support_point(n_components_solution);
    for (unsigned int j = 0; j < n_components_solution; ++j)
      solution_at_support_point[j] = solution[i * n_components_solution + j];

    // compute function at DoF support point
    std::vector<double> function_at_support_point =
      function(solution_at_support_point);
    for (unsigned int j = 0; j < n_components_function; ++j)
      function_dof_values[i * n_components_function + j] =
        function_at_support_point[j];
  }
}

/**
 * \brief Gets the function DoF values.
 *
 * \return function values at DoF support points
 */
template <int dim, bool is_scalar>
Vector<double> GroupFEValuesBase<dim, is_scalar>::get_function_dof_values() const
{
  // assert that function values are initialized
  Assert(function_values_initialized, ExcNotInitialized());

  return function_dof_values;
}
