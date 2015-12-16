/**
 * \file NoDiffusion.cc
 * \brief Provides the function definitions for the NoDiffusion class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
NoDiffusion<dim>::NoDiffusion()
  : ArtificialDiffusion<dim>()
{
}

/**
 * \brief Applies no viscosity (does nothing).
 *
 * \param[in] viscosity viscosity
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 * \param[in] fe_values FE values, reinitialized for cell already
 * \param[inout] cell_residual the residual vector for the cell
 */
/*
template <int dim>
void NoDiffusion<dim>::apply(std::shared_ptr<Viscosity<dim>>,
                             const Vector<double> &,
                             const Cell &,
                             const FEValues<dim> &,
                             Vector<double> &) const
{
}
*/

/**
 * \brief Computes artificial diffusion matrix.
 *
 * This returns a zero matrix.
 *
 * \param[in] viscosity pointer to viscosity cell map
 * \param[out] diffusion_matrix diffusion matrix
 */
template <int dim>
void NoDiffusion<dim>::compute_diffusion_matrix(
  const std::shared_ptr<Viscosity<dim>>, SparseMatrix<double> & diffusion_matrix)
{
  // return a zero matrix
  diffusion_matrix = 0;
}
