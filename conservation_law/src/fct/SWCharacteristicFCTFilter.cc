/**
 * \file SWCharacteristicFCTFilter.cc
 * \brief Provides the function definitions for the SWCharacteristicFCTFilter
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 * \param[in] gravity_  acceleration due to gravity \f$g\f$
 */
template <int dim>
SWCharacteristicFCTFilter<dim>::SWCharacteristicFCTFilter(
  const RunParameters & run_parameters_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const SparseMatrix<double> & lumped_mass_matrix_,
  const double & gravity_)
  : CharacteristicFCTFilter<dim>(
      run_parameters_, limiter_, dof_handler_, fe_, lumped_mass_matrix_),
    gravity(gravity_)
{
  // currently transformation matrix is limited to 1-D
  AssertThrow(dim == 1, ExcNotImplemented());
}

/**
 * \brief Computes a local characteristic transformation matrix.
 *
 * For the 1-D shallow water equations, the characteristic transformation matrix
 * is
 * \f[
 *   \mathbf{T}(\mathbf{u}) = \left[\begin{array}{c c}
 *     -\frac{2g}{a} & \frac{1}{h}\\
 *      \frac{2g}{a} & \frac{1}{h}
 *   \end{array}\right] \,,
 * \f]
 * where \f$u\f$ is the speed, and \f$a\f$ is the ``speed of sound''. This
 * matrix is defined such that
 * \f[
 *   \hat{\mathbf{u}} = \mathbf{T}(\mathbf{u})\mathbf{u} \,.
 * \f]
 *
 * \param[in] solution  solution vector \f$\mathbf{u}\f$ at which to evaluate
 *                      transformation
 *
 * \return transformation matrix \f$\mathbf{T}(\mathbf{u})\f$
 */
template <int dim>
FullMatrix<double> SWCharacteristicFCTFilter<dim>::compute_transformation_matrix(
  const Vector<double> & solution) const
{
  // extract solution components
  const double height = solution[0];

  // compute speed and sound speed
  const double a = std::sqrt(gravity * height);

  // compute matrix entries
  FullMatrix<double> matrix(this->n_components, this->n_components);
  matrix[0][0] = -2.0 * gravity / a;
  matrix[0][1] = 1.0 / height;
  matrix[1][0] = 2.0 * gravity / a;
  matrix[1][1] = 1.0 / height;

  return matrix;
}
