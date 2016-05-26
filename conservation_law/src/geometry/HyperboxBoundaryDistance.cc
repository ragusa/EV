/**
 * \file HyperboxBoundaryDistance.cc
 * \brief Provides the function definitions for the HyperboxBoundaryDistance
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] direction_  transport direction vector
 */
template <int dim>
HyperboxBoundaryDistance<dim>::HyperboxBoundaryDistance(
  const Tensor<1, dim> & direction_,
  const double & x_min_,
  const double & y_min_,
  const double & z_min_)
  : BoundaryDistance<dim>(direction_)
{
  // initialize hyperbox lower boundary in each dimension
  x_min.resize(dim);
  x_min[0] = x_min_;
  if (dim > 1)
    x_min[1] = y_min_;
  if (dim > 2)
    x_min[2] = z_min_;
}

/**
 * \brief Computes distance to boundary from a point.
 *
 * For a hyperbox (aligned with axes) with lower boundaries
 * \f$(x^-_1,x^-_2,x^-_3)\f$, the distance is computed as
 * \f[
 *   s = \sqrt{s_i^2 + s_j^2 + s_k^2} ,
 * \f]
 * where \f$i\f$ corresponds to the index for which the following is true:
 * \f[
 *   i \quad : \quad \frac{x_i - x^-_i}{\Omega_i}
 *     = \min\limits_j\left(\frac{x_j - x^-_j}{\Omega_j}\right) ,
 * \f]
 * where \f$\mathbf{\Omega} \equiv (\Omega_1,\Omega_2,\Omega_3)\f$ is the
 * direction vector, and \f$s_j\f$ are computed as
 * \f[
 *   s_i = x_i - x_i^- , \quad
 *   s_j = \frac{\Omega_j}{\Omega_i}s_i , \quad
 *   s_k = \frac{\Omega_k}{\Omega_i}s_i .
 * \f]
 *
 * \param[in] p  point at which to evaluate distance to boundary
 *
 * return exact solution value at point
 */
template <int dim>
double HyperboxBoundaryDistance<dim>::compute(const Point<dim> & p) const
{
  // components of distance travelled in region
  std::vector<double> s(3, 0.0);

  // index for direction that has the minimum distance ratio
  unsigned int d_min = 0;

  // initialize to x-direction
  s[0] = p[0] - x_min[0];
  double r_min = s[0] / this->direction[0];

  // loop over the number of dimensions and compare distance ratios
  // to determine which absorber region plane segment through which the
  // particle passes
  for (unsigned int d = 1; d < dim; ++d)
  {
    // compute distance ratio for this direction
    double r_d = (p[d] - x_min[d]) / this->direction[d];

    // compare to the minimum distance comp
    if (r_d < r_min)
    { // this direction is new minimum distance ratio

      // set minimum distance ratio index to this direction
      d_min = d;

      // set minimum distance ratio to this one
      r_min = r_d;

      // compute distance for this direction
      s[d] = p[d] - x_min[d];

      // recompute distance for all previous directions
      for (unsigned int dd = 0; dd <= d - 1; ++dd)
        s[dd] = this->direction[dd] / this->direction[d] * s[d];
    }
    else
    {
      // compute distance for this direction
      s[d] = this->direction[d] / this->direction[d_min] * s[d_min];
    }
  }

  // compute total distance
  double total_distance;
  if (dim == 1)
  {
    total_distance = s[0] / this->direction[0];
  }
  else if (dim == 2)
  {
    // compute direction component in z direction using the fact that direction
    // vector is normalized to 1: omega_x^2 + omega_y^2 + omega_z^2 = 1
    const double direction_z2 =
      1.0 - std::pow(this->direction[0], 2) - std::pow(this->direction[1], 2);
    const double sz2 =
      direction_z2 * std::pow(s[d_min] / this->direction[d_min], 2);
    total_distance = std::sqrt(std::pow(s[0], 2) + std::pow(s[1], 2) + sz2);
  }
  else // dim == 3
  {
    total_distance =
      std::sqrt(std::pow(s[0], 2) + std::pow(s[1], 2) + std::pow(s[2], 2));
  }

  return total_distance;
}
