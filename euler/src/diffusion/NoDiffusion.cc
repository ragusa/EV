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
 * \brief Computes artificial diffusion matrix.
 *
 * This returns a zero matrix.
 *
 * \param[in] solution solution vector
 * \param[in] viscosity pointer to viscosity cell map
 * \param[out] diffusion_matrix diffusion matrix
 */
template <int dim>
void NoDiffusion<dim>::compute_diffusion_matrix(
  const Vector<double> &,
  const std::shared_ptr<Viscosity<dim>>,
  SparseMatrix<double> & diffusion_matrix)
{
  // return a zero matrix
  diffusion_matrix = 0;
}
