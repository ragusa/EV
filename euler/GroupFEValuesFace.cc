/**
 * \file GroupFEValuesFace.cc
 * \brief Provides the function definitions for the GroupFEValuesFace class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] n_components_solution_ number of solution components
 * \param[in] n_components_function_ number of function components
 * \param[in] solution_dof_handler_ DoF handler for solution
 * \param[in] triangulation_ triangulation
 * \param[in] face_quadrature_ face quadrature
 * \param[in] solution_ solution vector
 * \param[in] aux_vector_ optional auxiliary vector
 */
template <int dim, bool is_scalar>
GroupFEValuesFace<dim, is_scalar>::GroupFEValuesFace(
  const unsigned int & n_components_solution_,
  const unsigned int & n_components_function_,
  const DoFHandler<dim> & solution_dof_handler_,
  const Triangulation<dim> & triangulation_,
  const QGauss<dim-1> & face_quadrature_,
  const Vector<double> & solution_,
  const Vector<double> & aux_vector_)
  : GroupFEValuesBase<dim, is_scalar>(n_components_solution_,
                                      n_components_function_,
                                      solution_dof_handler_,
                                      triangulation_,
                                      solution_,
                                      aux_vector_),
    face_quadrature(face_quadrature_),
    n_quadrature_points(face_quadrature.size()),
    function_fe_values_face(this->fe,
                            face_quadrature,
                            update_values | update_gradients |
                              update_normal_vectors)
{
}

/**
 * \brief Reinitializes function FE values for face.
 *
 * \param[in] solution_cell solution cell iterator
 * \param[in] face face index
 */
template <int dim, bool is_scalar>
void GroupFEValuesFace<dim, is_scalar>::reinit(const Cell & solution_cell,
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
template <int dim, bool is_scalar>
void GroupFEValuesFace<dim, is_scalar>::get_function_values(
  std::vector<ValueType> & function_values) const
{
  function_fe_values_face[this->function_extractor].get_function_values(
    this->function_dof_values, function_values);
}

/**
 * \brief Computes function gradients at quadrature points in a face.
 *
 * \param[out] function_gradients vector of function gradients at quadrature
 *             points in a face
 */
template <int dim, bool is_scalar>
void GroupFEValuesFace<dim, is_scalar>::get_function_gradients(
  std::vector<GradientType> & function_gradients) const
{
  function_fe_values_face[this->function_extractor].get_function_gradients(
    this->function_dof_values, function_gradients);
}

/**
 * \brief Gets normal vectors at quadrature points in a face.
 *
 * \return normal vectors at quadrature points in a face
 */
template <int dim, bool is_scalar>
std::vector<Tensor<1,dim>> GroupFEValuesFace<dim, is_scalar>::get_normal_vectors() const
{
  return function_fe_values_face.get_all_normal_vectors();
}
