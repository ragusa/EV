/**
 * \file GroupFEValuesFace.cc
 * \brief Provides the function definitions for the GroupFEValuesFace class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_components_ number of solution components
 * \param[in] solution_dof_handler_ DoF handler for solution
 * \param[in] triangulation_ triangulation
 * \param[in] face_quadrature_ face quadrature
 * \param[in] solution_ solution vector
 * \param[in] aux_vector_ optional auxiliary vector
 */
template <int dim>
GroupFEValuesFace<dim>::GroupFEValuesFace(
  const unsigned int & n_components_,
  const DoFHandler<dim> & solution_dof_handler_,
  const Triangulation<dim> & triangulation_,
  const QGauss<dim> & face_quadrature_,
  const Vector<double> & solution_,
  const Vector<double> & aux_vector_)
  : GroupFEValuesBase<dim>(
      n_components_, solution_dof_handler_, triangulation_, solution_, aux_vector_),
    face_quadrature(face_quadrature_),
    n_quadrature_points(face_quadrature.size()),
    function_fe_values_face(
      this->fe, face_quadrature, update_values | update_gradients)
{
}

/**
 * \brief Reinitializes function FE values for face.
 *
 * \param[in] solution_cell solution cell iterator
 * \param[in] face face index
 */
template <int dim>
void GroupFEValuesFace<dim>::reinit(const cell_iterator & solution_cell,
                                    const unsigned int & face)
{
  // reinitialize function FE values with the function cell corresponding
  // to the solution cell
  function_fe_values_face.reinit(
    this->solution_cell_to_function_cell_map[solution_cell], face);
}

/**
 * \brief Computes function values at quadrature points in a face.
 *
 * \param[out] function_values vector of function values at quadrature points
 *             in a face
 */
template <int dim>
void GroupFEValuesFace<dim>::get_function_values(
  std::vector<double> & function_values) const
{
  function_fe_values_face.get_function_values(this->function_dof_values,
                                              function_values);
}
