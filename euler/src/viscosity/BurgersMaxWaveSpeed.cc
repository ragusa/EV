/**
 * \file BurgersMaxWaveSpeed.cc
 * \brief Provides the function definitions for the BurgersMaxWaveSpeed
 *        class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
BurgersMaxWaveSpeed<dim>::BurgersMaxWaveSpeed()
  : MaxWaveSpeed<dim>()
{
}

/**
 * \brief Computes the maximum wave speed for a given left and right state.
 *
 * \param[in] solution_left solution vector for left state
 * \param[in] solution_right solution vector for right state
 * \param[in] normal_vector normal vector
 *
 * \return maximum wave speed
 */
template <int dim>
double BurgersMaxWaveSpeed<dim>::compute(
  const std::vector<double> & solution_left,
  const std::vector<double> & solution_right,
  const Tensor<1, dim> &) const
{
  // extract solution components from solution vector
  const double uL = solution_left[0];
  const double uR = solution_right[0];

  // compute maximum wave speed
  double max_wave_speed;
  if (uL <= uR) // rarefaction
    max_wave_speed = std::max(std::abs(uL), std::abs(uR));
  else // shock
    max_wave_speed = 0.5 * (uR * uR - uL * uL) / (uR - uL);

  return max_wave_speed;
}
