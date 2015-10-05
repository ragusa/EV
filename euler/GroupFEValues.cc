/**
 * \file GroupFEValues.cc
 * \brief Provides the function definitions for the GroupFEValues class.
 */

/**
 * \brief Constructor.
 *
 * In addition to construction, this function distributes DoFs and initializes
 * function cell iterators. Derived classes must compute the function
 * DoF values (it is convenient to do this in their constructors).
 *
 * \param[in] n_components number of solution components
 * \param[in] solution solution vector
 * \param[in] triangulation triangulation
 * \param[in] cell_quadrature cell quadrature
 */
template <int dim>
GroupFEValues<dim>::GroupFEValues(
  const unsigned int & n_components_,
  const Vector<double> & solution_,
  const Triangulation<dim> & triangulation_,
  const QGauss<dim> & cell_quadrature_) :
  n_components(n_components_),
  fe(1),
  dof_handler(triangulation_),
  n_function_dofs_per_cell(fe.dofs_per_cell),
  n_solution_dofs_per_cell(n_function_dofs_per_cell*n_components),
  cell_quadrature(cell_quadrature_),
  n_quadrature_points(cell_quadrature.size()),
  solution(&solution_),
  fe_values_function(fe, cell_quadrature, update_values | update_gradients)
{
  // distribute DoFs and initialize cell iterators
  dof_handler.distribute_dofs(fe);
  n_dofs = dof_handler.n_dofs();
  cell_function = dof_handler.begin_active();
  endc_function = dof_handler.end();

  // resize function DoFs vector
  function_dof_values.resize(n_dofs);
}

/**
 * \brief Computes function values at quadrature points and advances
 *        the function cell iterator.
 */
template <int dim>
std::vector<double> GroupFEValues<dim>::get_function_values(
  const typename DoFHandler<dim>::active_cell_iterator & cell)
{
  // get local solution DoF indices and determine local function DoF indices
  // and then the corresponding function DoFs
  std::vector<unsigned int> local_dof_indices(n_solution_dofs_per_cell);
  cell->get_dof_indices(local_dof_indices);
  std::vector<double> function_dofs_local(n_function_dofs_per_cell);
  for (unsigned int i = 0; i < n_function_dofs_per_cell; ++i)
  {
    const unsigned int k = local_dof_indices[i*n_components];
    Assert(k % n_components == 0, ExcModulusNotZero(k,n_components));
    const unsigned int function_dof_index = k / n_components;
    function_dofs_local[i] = function_dof_values[function_dof_index];
  }

  // reinitialize finite element values for function
  fe_values_function.reinit(cell_function);

  // compute function values at the quadrature points
  std::vector<double> function_values(n_quadrature_points, 0.0);
  for (unsigned int q = 0; q < n_quadrature_points; ++q)
    for (unsigned int i = 0; i < n_function_dofs_per_cell; ++i)
      function_values[q] +=
        function_dofs_local[i]*fe_values_function.shape_value(i, q);

  // advance function cell iterator
  cell_function++;

  return function_values;
}

/**
 * \brief Computes the function DoF values.
 */
template <int dim>
void GroupFEValues<dim>::compute_function_dof_values()
{
  function_dof_values.resize(n_dofs);
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // extract all components of solution at DoF support point
    std::vector<double> solution_at_support_point(n_components);
    for (unsigned int j = 0; j < n_components; ++j)
      solution_at_support_point[j] = (*solution)[i*n_components + j];

    // compute function at DoF support point
    function_dof_values[i] = function(solution_at_support_point);
  }
}

/**
 * \brief Gets the function DoF values.
 *
 * \return function values at DoF support points
 */
template <int dim>
std::vector<double> GroupFEValues<dim>::get_function_dof_values() const
{
  return function_dof_values;
}
