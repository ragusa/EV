/**
 * \file ExplicitEulerFCT.cc
 * \brief Provides the function definitions for the ExplicitEulerFCT class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ExplicitEulerFCT<dim>::ExplicitEulerFCT()
{
}

/**
 * \brief Computes the limited antidiffusion vector with entries
 * \f$\bar{p}_i=\sum\limits_j L_{i,j}P_{i,j}\f$.
 *
 * \param[out] antidiffusion_vector  the antidiffusion vector
 *             \f$\bar{\mathbf{p}}\f$
 */
template <int dim>
void ExplicitEulerFCT<dim>::compute_antidiffusion_vector_fe(
  Vector<double> & antidiffusion_vector)
{
  // compute antidiffusive fluxes
  compute_antidiffusive_fluxes_fe();

  // filter the antidiffusive fluxes
  this->filter_antidiffusive_fluxes();

  // compute sum of limited antidiffusive fluxes
  this->compute_limited_antidiffusion_sum(antidiffusion_vector);
}

