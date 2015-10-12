/**
 * \file ShallowWaterRiemannSolver.cc
 * \brief Provides the function definitions for the ShallowWaterRiemannSolver
 *        class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ShallowWaterRiemannSolver<dim>::ShallowWaterRiemannSolver(
  const double & h_left_,
  const double & u_left_,
  const double & h_right_,
  const double & u_right_,
  const double & gravity_,
  const double & x_interface_)
  : Function<dim>(dim + 1),
    h_left(h_left_),
    u_left(u_left_),
    h_right(h_right_),
    u_right(u_right_),
    gravity(gravity_),
    x_interface(x_interface_),
    c_left(std::sqrt(gravity * h_left)),
    c_right(std::sqrt(gravity * h_right))
{
  // assert that this is 1-D
  Assert(dim == 1, ExcImpossibleInDim(dim));

  // compute initial guess for height - start with 2-rarefaction solution
  h_star =
    std::pow(0.5 * (c_left + c_right) - 0.25 * (u_right - u_left), 2) / gravity;

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
    evaluate_fluxes(h_star_old, h_left, c_left, f_left, f_deriv_left);
    evaluate_fluxes(h_star_old, h_right, c_right, f_right, f_deriv_right);

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

  // compute velocity and sound speed in star region
  u_star = 0.5 * (u_left + u_right) + 0.5 * (f_right - f_left);
  c_star = std::sqrt(gravity * h_star);
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
void ShallowWaterRiemannSolver<dim>::evaluate_fluxes(const double & h,
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

/**
 * \brief Computes exact solution at a single point.
 *
 * \param[in] p point at which to evaluate exact solution
 * \param[in] component component index for which to evaluate
 *
 * \return solution for a component at a point
 */
template <int dim>
double ShallowWaterRiemannSolver<dim>::value(const Point<dim> & p,
                                             const unsigned int component) const
{
  // get time
  double t = this->get_time();

  // get x position
  double x = p[0];

  // sample the solution
  double h_sample;
  double u_sample;
  if (t == 0)
  {
    if (x <= x_interface)
    {
      h_sample = h_left;
      u_sample = u_left;
    }
    else
    {
      h_sample = h_right;
      u_sample = u_right;
    }
  }
  else
  {
    // sample point
    double s = (x - x_interface) / t;

    if (s <= u_star) // left of contact discontinuity
    {
      if (h_star >= h_left) // left shock
      {
        // compute left wave speed
        double q_left =
          std::sqrt(((h_star + h_left) * h_star) / (2.0 * h_left * h_left));
        double s_left = u_left - c_left * q_left;

        if (s <= s_left) // left of shock; left state
        {
          h_sample = h_left;
          u_sample = u_left;
        }
        else // right of shock; star state
        {
          h_sample = h_star;
          u_sample = u_star;
        }
      }
      else // left rarefaction
      {
        // compute wave speed of head of rarefaction
        double s_head_left = u_left - c_left;

        if (s <= s_head_left) // left of rarefaction; left state
        {
          h_sample = h_left;
          u_sample = u_left;
        }
        else
        {
          // compute wave speed of tail of rarefaction
          double s_tail_left = u_star - c_star;

          if (s <= s_tail_left) // inside rarefaction
          {
            u_sample = (u_left + 2.0 * c_left + 2.0 * s) / 3.0;
            double c_sample = (u_left + 2.0 * c_left - s) / 3.0;
            h_sample = std::pow(c_sample, 2) / gravity;
          }
          else // right of rarefaction; star state
          {
            h_sample = h_star;
            u_sample = u_star;
          }
        }
      }
    }
    else // right of contact discontinuity
    {
      if (h_star >= h_right) // right shock
      {
        // compute right wave speed
        double q_right =
          std::sqrt(((h_star + h_right) * h_star) / (2.0 * h_right * h_right));
        double s_right = u_right + c_right * q_right;

        if (s >= s_right) // right of shock; right state
        {
          h_sample = h_right;
          u_sample = u_right;
        }
        else // left of shock; star state
        {
          h_sample = h_star;
          u_sample = u_star;
        }
      }
      else // right rarefaction
      {
        // compute wave speed of head of rarefaction
        double s_head_right = u_right + c_right;

        if (s >= s_head_right) // right of rarefaction; right state
        {
          h_sample = h_right;
          u_sample = u_right;
        }
        else
        {
          // compute wave speed of tail of rarefaction
          double s_tail_right = u_star + c_star;

          if (s >= s_tail_right) // inside rarefaction
          {
            u_sample = (u_right - 2.0 * c_right + 2.0 * s) / 3.0;
            double c_sample = (-u_right + 2.0 * c_right + s) / 3.0;
            h_sample = std::pow(c_sample, 2) / gravity;
          }
          else // left of rarefaction; star state
          {
            h_sample = h_star;
            u_sample = u_star;
          }
        }
      }
    }
  }

  // return solution for the requested component
  if (component == 0)
    return h_sample;
  else
    return h_sample * u_sample;
}
