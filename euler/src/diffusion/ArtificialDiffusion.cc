/**
 * \file ArtificialDiffusion.cc
 * \brief Provides the function definitions for the ArtificialDiffusion class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ArtificialDiffusion<dim>::ArtificialDiffusion()
{
}

/**
 * \brief Reinitializes diffusion matrix.
 *
 * \param[in] sparsity sparsity pattern
 * \param[in] dof_handler degree of freedom handler
 */
template <int dim>
void ArtificialDiffusion<dim>::reinitialize(const SparsityPattern & sparsity,
                                            const DoFHandler<dim> &)
{
  diffusion_matrix.reinit(sparsity);
}

/**
 * \brief Updates the diffusion matrix.
 *
 * Here the function does nothing, but derived classes can override it.
 *
 * \param[in] solution solution vector
 */
template <int dim>
void ArtificialDiffusion<dim>::update(const Vector<double> & solution)
{
}

/**
 * \brief Applies diffusion matrix product with solution if the artificial
 *        diffusion is of the algebraic type.
 *
 * Here the function does nothing, but derived classes can override it.
 *
 * \param[in] solution solution vector
 * \param[inout] ss_residual steady-state residual
 */
template <int dim>
void ArtificialDiffusion<dim>::apply_algebraic_diffusion(const Vector<double> &,
                                                         Vector<double> &)
{
}
