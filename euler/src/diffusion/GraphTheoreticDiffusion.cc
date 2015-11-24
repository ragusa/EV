/**
 * \file GraphTheoreticDiffusion.cc
 * \brief Provides the function definitions for the GraphTheoreticDiffusion class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] dofs_per_cell_ number of degrees of freedom per cell
 * \param[in] n_components_ number of solution components
 */
template <int dim>
GraphTheoreticDiffusion<dim>::GraphTheoreticDiffusion(
  const unsigned int & dofs_per_cell_, const unsigned int & n_components_)
  : ArtificialDiffusion<dim>(),
    dofs_per_cell(dofs_per_cell_),
    n_components(n_components_),
    dofs_per_cell_per_component(dofs_per_cell_ / n_components_)
{
}

/**
 * \brief Applies graph-theoretic diffusion.
 *
 * \param[in] viscosity viscosity
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 * \param[in] fe_values FE values, reinitialized for cell already
 * \param[inout] cell_residual the residual vector for the cell
 */
template <int dim>
void GraphTheoreticDiffusion<dim>::apply(
  std::shared_ptr<Viscosity<dim>> viscosity,
  const Vector<double> & solution,
  const Cell & cell,
  const FEValues<dim> & fe_values,
  Vector<double> & cell_residual) const;
{
  // get cell volume
  double cell_volume = cell->measure();

  // loop over DoFs in cell
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
  {
    // compute b_K(u,test[i])
    double b_i = 0.0;
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      // if j corresponds to the same solution component as i
      // The ordering is: u^0[0],u^1[0],u^0[1],u^1[1],etc.
      // therefore DoF indices of the same component differ by a factor
      // of n_components
      if ((i - j) % n_components == 0) // if m(i) == m(j)
      {
        // compute local viscous bilinear form: b_K(test[j],test[i])
        double b_K;
        if (j == i)
          b_K = cell_volume;
        else
          b_K = -1.0 / (dofs_per_cell_per_component - 1.0) * cell_volume;

        // add U[j]*b_K(test[j],test[i]) to sum for b_K(u,test[i])
        b_i += solution(local_dof_indices[j]) * b_K;
      }
    }
    cell_residual(i) -= viscosity[cell] * b_i;
  }
}
