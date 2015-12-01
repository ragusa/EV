/**
 * \file ShallowWaterMaxWaveSpeed.cc
 * \brief Provides the function definitions for the ShallowWaterMaxWaveSpeed
 *        class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] gravity_ acceleration due to gravity
 */
template <int dim>
ShallowWaterMaxWaveSpeed<dim>::ShallowWaterMaxWaveSpeed(const double & gravity_)
  : MaxWaveSpeed<dim>(), gravity(gravity_)
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
double ShallowWaterMaxWaveSpeed<dim>::compute(
  const std::vector<double> & solution_left,
  const std::vector<double> & solution_right,
  const Tensor<1, dim> & normal_vector) const
{
  // extract solution components from solution vector
  const double height_left = solution_left[0];
  const double height_right = solution_right[0];
  Tensor<1, dim> momentum_left_multidim;
  Tensor<1, dim> momentum_right_multidim;
  for (unsigned int d = 0; d < dim; ++d)
  {
    momentum_left_multidim[d] = solution_left[d + 1];
    momentum_right_multidim[d] = solution_right[d + 1];
  }

  // compute normal component of velocities
  const double velocity_left =
    (momentum_left_multidim * normal_vector) / height_left;
  const double velocity_right =
    (momentum_right_multidim * normal_vector) / height_right;

  // compute height in star region
  const double height_star =
    compute_height_star(height_left, velocity_left, height_right, velocity_right);

  // compute sound speeds
  const double soundspeed_left = std::sqrt(gravity * height_left);
  const double soundspeed_right = std::sqrt(gravity * height_right);

  // compute left and right head speeds
  const double left_head_speed = velocity_left -
    soundspeed_left * (1.0 + std::max(0.5 * (height_star - height_left) *
                                        (height_star + 2.0 * height_left) /
                                        std::pow(height_left, 2),
                                      0.0));
  const double right_head_speed = velocity_right +
    soundspeed_right * (1.0 + std::max(0.5 * (height_star - height_right) *
                                         (height_star + 2.0 * height_right) /
                                         std::pow(height_right, 2),
                                       0.0));

  // compute maximum wave speed
  const double max_wave_speed = std::max(left_head_speed, right_head_speed);

  return max_wave_speed;
}

/**
 * \brief Computes height in star region.
 *
 * \param[in] h_left height in left state
 * \param[in] u_left velocity in left state
 * \param[in] h_right height in right state
 * \param[in] u_right velocity in right state
 *
 * \return height in star region
 */
template <int dim>
double ShallowWaterMaxWaveSpeed<dim>::compute_height_star(
  const double & h_left,
  const double & u_left,
  const double & h_right,
  const double & u_right) const
{
  // compute sound speeds
  const double a_left = std::sqrt(gravity * h_left);
  const double a_right = std::sqrt(gravity * h_right);

  // compute initial guess for height - start with 2-rarefaction solution
  double h_star =
    std::pow(0.5 * (a_left + a_right) - 0.25 * (u_right - u_left), 2) / gravity;

  double h_min = std::min(h_left, h_right);
  if (h_star <= h_min)
  {
    // keep 2-rarefaction solution
  }
  else
  {
    // use 2-shock solution
    double GE_left =
      std::sqrt(0.5 * gravity * (h_star + h_left) / (h_star * h_left));
    double GE_right =
      std::sqrt(0.5 * gravity * (h_star + h_right) / (h_star * h_right));
    h_star = (GE_left * h_left + GE_right * h_right - (u_right - u_left)) /
      (GE_left + GE_right);
  }

  // begin Newton-Raphson iteration to determine h in star region
  const unsigned int max_iteration = 50; // maximum number of iterations
  const double tolerance = 1.0e-6;       // iteration tolerance
  double h_star_old = h_star;            // old iterate for h in star region
  bool converged;                        // flag for convergence
  double f_left = 0.0, f_right = 0.0, f_deriv_left = 0.0,
         f_deriv_right = 0.0; // fluxes
  for (unsigned int k = 0; k < max_iteration; ++k)
  {
    // evaluate left and right fluxes and their derivatives
    evaluate_fluxes(h_star_old, h_left, a_left, f_left, f_deriv_left);
    evaluate_fluxes(h_star_old, h_right, a_right, f_right, f_deriv_right);

    // compute new iterate for h in star region
    h_star = h_star_old -
      (f_left + f_right + u_right - u_left) / (f_deriv_left + f_deriv_right);

    // compute relative change
    double dh = std::abs(h_star - h_star_old) / (0.5 * (h_star + h_star_old));

    // determine if iteration has converged
    if (dh <= tolerance)
    {
      converged = true;
      break;
    }
    else
    {
      converged = false;
    }

    // cancel negative h
    if (h_star <= 0.0)
    {
      h_star = tolerance;
    }

    // reset for next iteration
    h_star_old = h_star;
  }

  // check to ensure iteration converged
  Assert(converged, ExcInvalidState());

  return h_star;
}

/**
 * \brief Evaluates flux and its derivative.
 *
 * \param[in] h height at which to evaluate flux
 * \param[in] h_K height of side K
 * \param[in] c_K sound speed of side K
 * \param[out] f flux
 * \param[out] f_deriv flux derivative
 */
template <int dim>
void ShallowWaterMaxWaveSpeed<dim>::evaluate_fluxes(const double & h,
                                                    const double & h_K,
                                                    const double & c_K,
                                                    double & f,
                                                    double & f_deriv) const
{
  if (h <= h_K)
  {
    // wave is rarefaction
    double c = std::sqrt(gravity * h);
    f = 2.0 * (c - c_K);
    f_deriv = gravity / c;
  }
  else
  {
    // wave is shock
    double GE_star = std::sqrt(0.5 * gravity * (h + h_K) / (h * h_K));
    f = (h - h_K) * GE_star;
    f_deriv = GE_star - (0.25 * gravity * (h - h_K)) / (GE_star * h * h);
  }
}
