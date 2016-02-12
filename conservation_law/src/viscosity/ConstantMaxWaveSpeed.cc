/**
 * \file ConstantMaxWaveSpeed.cc
 * \brief Provides the function definitions for the ConstantMaxWaveSpeed
 *        class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] constant_speed_ constant maximum wave speed
 */
template <int dim>
ConstantMaxWaveSpeed<dim>::ConstantMaxWaveSpeed(const double & constant_speed_)
  : MaxWaveSpeed<dim>(), constant_speed(constant_speed_)
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
double ConstantMaxWaveSpeed<dim>::compute(const std::vector<double> &,
                                          const std::vector<double> &,
                                          const Tensor<1, dim> &) const
{
  return constant_speed;
}
