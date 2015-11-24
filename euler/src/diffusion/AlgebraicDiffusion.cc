/**
 * \file AlgebraicDiffusion.cc
 * \brief Provides the function definitions for the AlgebraicDiffusion class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
AlgebraicDiffusion<dim>::AlgebraicDiffusion()
{
}

/**
 * \brief Reinitializes diffusion matrix and temporary vector.
 *
 * \param[in] sparsity sparsity pattern
 * \param[in] dof_handler degree of freedom handler
 */
template <int dim>
void AlgebraicDiffusion<dim>::reinitialize(const SparsityPattern & sparsity,
                                           const DoFHandler<dim> & dof_handler)
{
  // reinitialize diffusion matrix
  diffusion_matrix.reinit(sparsity);

  // resize temporary vector
  tmp_vector.reinit(dof_handler.n_dofs());
}

/**
 * \brief Does nothing; application is performed after cell loop.
 *
 * \param[in] viscosity viscosity
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 * \param[in] fe_values FE values, reinitialized for cell already
 * \param[inout] cell_residual the residual vector for the cell
 */
template <int dim>
void AlgebraicDiffusion<dim>::apply(std::shared_ptr<Viscosity<dim>>,
                                    const Vector<double> &,
                                    const Cell &,
                                    const FEValues<dim> &,
                                    Vector<double> &) const
{
}

/**
 * \brief Applies diffusion matrix product with solution.
 *
 * \param[in] solution solution vector
 * \param[inout] ss_residual steady-state residual
 */
template <int dim>
void AlgebraicDiffusion<dim>::apply_algebraic_diffusion(
  const Vector<double> & solution, Vector<double> & ss_residual)
{
  // compute product of diffusion matrix and solution
  this->diffusion_matrix.vmult(tmp_vector, solution);

  // subtract result from steady-state residual
  ss_residual.add(-1.0, tmp_vector);
}
