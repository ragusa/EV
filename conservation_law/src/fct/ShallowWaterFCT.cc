/**
 * \file ShallowWaterFCT.cc
 * \brief Provides the function definitions for the ShallowWaterFCT class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] parameters_  parameters
 * \param[in] dof_handler_  DoF handler
 * \param[in] triangulation_  triangulation
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 * \param[in] consistent_mass_matrix_  consistent mass matrix \f$\mathbf{M}^C\f$
 * \param[in] star_state_  star state object for computing
 *                         \f$\mathbf{U}^*_{i,j}\f$
 * \param[in] linear_solver_  linear solver
 * \param[in] sparsity_pattern_  sparsity pattern
 * \param[in] dirichlet_nodes_  vector of Dirichlet nodes
 * \param[in] n_components_  number of components
 * \param[in] dofs_per_cell_  number of degrees of freedom per cell
 * \param[in] component_names_  list of component names
 * \param[in] use_star_states_in_fct_bounds_  flag to signal that star states are
 *            to be used in FCT bounds
 * \param[in] gravity_  acceleration due to gravity
 */
template <int dim>
ShallowWaterFCT<dim>::ShallowWaterFCT(
  const RunParameters<dim> & parameters_,
  const DoFHandler<dim> & dof_handler_,
  const Triangulation<dim> & triangulation_,
  const SparseMatrix<double> & lumped_mass_matrix_,
  const SparseMatrix<double> & consistent_mass_matrix_,
  const std::shared_ptr<StarState<dim>> & star_state_,
  const LinearSolver<dim> & linear_solver_,
  const SparsityPattern & sparsity_pattern_,
  const std::vector<unsigned int> & dirichlet_nodes_,
  const unsigned int & n_components_,
  const unsigned int & dofs_per_cell_,
  const std::vector<std::string> & component_names_,
  const bool & use_star_states_in_fct_bounds_,
  const double & gravity_)
  : FCT<dim>(parameters_,
             dof_handler_,
             triangulation_,
             lumped_mass_matrix_,
             consistent_mass_matrix_,
             star_state_,
             linear_solver_,
             sparsity_pattern_,
             dirichlet_nodes_,
             n_components_,
             dofs_per_cell_,
             component_names_,
             use_star_states_in_fct_bounds_),
    gravity(gravity_)
{
}

/**
 * \brief Computes a local characteristic transformation matrix.
 *
 * For the shallow water equations, the characteristic transformation matrix is
 * \f[
 *   \mathbf{T}(\mathbf{u}) = \left[\begin{array}{c c}
 *     1   & 1\\
 *     u-a & u+a
 *   \end{array}\right] \,,
 * \f]
 * where \f$u\f$ is the speed, and \f$a\f$ is the ``speed of sound''
 *
 * \param[in] solution  solution vector \f$\mathbf{u}\f$ at which to evaluate
 *                      transformation
 *
 * \return transformation matrix \f$\mathbf{T}(\mathbf{u})\f$
 */
template <int dim>
FullMatrix<double> ShallowWaterFCT<dim>::compute_transformation_matrix(
  const Vector<double> & solution) const
{
  Assert(dim == 1, ExcNotImplemented());

  // extract solution components
  const double height = solution[0];
  std::vector<double> momentum(1);
  for (unsigned int d = 0; d < dim; ++d)
    momentum[d] = solution[d + 1];

  // compute speed and sound speed
  const double u = momentum[0] / height;
  const double a = std::sqrt(gravity * height);

  // compute matrix entries
  FullMatrix<double> matrix(this->n_components, this->n_components);
  matrix[0][0] = 1.0;
  matrix[0][1] = 1.0;
  matrix[1][0] = u - a;
  matrix[1][1] = u + a;

  return matrix;
}

/**
 * \brief Computes a local characteristic transformation matrix inverse.
 *
 * For the shallow water equations, the characteristic transformation matrix
 * inverse is
 * \f[
 *   \mathbf{T}^{-1}(\mathbf{u}) = \frac{1}{2a}\left[\begin{array}{c c}
 *     u+a & -1\\
 *     a-u & 1
 *   \end{array}\right] \,,
 * \f]
 * where \f$u\f$ is the speed, and \f$a\f$ is the ``speed of sound''
 *
 * \param[in] solution  solution vector \f$\mathbf{u}\f$ at which to evaluate
 *                      transformation
 *
 * \return transformation matrix inverse \f$\mathbf{T}^{-1}(\mathbf{u})\f$
 */
template <int dim>
FullMatrix<double> ShallowWaterFCT<dim>::compute_transformation_matrix_inverse(
  const Vector<double> & solution) const
{
  Assert(dim == 1, ExcNotImplemented());

  // extract solution components
  const double height = solution[0];
  std::vector<double> momentum(1);
  for (unsigned int d = 0; d < dim; ++d)
    momentum[d] = solution[d + 1];

  // compute speed and sound speed
  const double u = momentum[0] / height;
  const double a = std::sqrt(gravity * height);

  // compute matrix entries
  FullMatrix<double> matrix(this->n_components, this->n_components);
  const double factor = 0.5 / a;
  matrix[0][0] = (u + a) * factor;
  matrix[0][1] = -factor;
  matrix[1][0] = (a - u) * factor;
  matrix[1][1] = factor;

  return matrix;
}
