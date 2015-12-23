/**
 * \file MaxWaveSpeed.cc
 * \brief Provides the function definitions for the MaxWaveSpeed class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
MaxWaveSpeed<dim>::MaxWaveSpeed(
  const std::shared_ptr<StarState<dim>> & star_state_)
  : max_wave_speed_domain(0.0), star_state(star_state_)
{
}

/**
 * \brief Resets the maximum wave speed in domain to zero.
 *
 * This function should be called at the beginning of each time step.
 */
template <int dim>
void MaxWaveSpeed<dim>::reset_max_wave_speed_domain()
{
  max_wave_speed_domain = 0.0;
}

/**
 * \brief Returns the maximum wave speed in the domain.
 *
 * \return maximum wave speed in the domain.
 */
template <int dim>
double MaxWaveSpeed<dim>::get_max_wave_speed_domain() const
{
  return max_wave_speed_domain;
}
