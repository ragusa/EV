/**
 * \file Entropy.cc
 * \brief Provides the function definitions for the Entropy class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
Entropy<dim>::Entropy()
  : need_to_compute_average_entropy()
{
}

/**
 * \brief Reinitializes group FE values (if any) and computes
 *        solution-dependent quantities for normalization of entropy
 *        viscosity.
 *
 * \param[in] solution solution vector
 */
void Entropy<dim>::reinitialize(const Vector<double> & solution)
{
  // reinitialize group FE values (if any)
  reinitialize_group_fe_values(solution);

  // compute average entropy if needed for normalization of entropy viscosity
  if (need_to_compute_average_entropy)
    entropy_average = compute_average_entropy(solution);
}
