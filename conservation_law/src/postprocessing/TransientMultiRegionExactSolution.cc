/**
 * \brief Constructor.
 *
 * \param[in] interface_positions_  positions interfaces between regions
 * \param[in] source_  source values in each region
 * \param[in] sigma_  sigma values in each region
 * \param[in] direction_  transport direction vector
 * \param[in] incoming_  incoming boundary value
 */
template <int dim>
TransientMultiRegionExactSolution<dim>::TransientMultiRegionExactSolution(
  const std::vector<double> & interface_positions_,
  const std::vector<double> & source_,
  const std::vector<double> & sigma_,
  const Tensor<1, dim> & direction_,
  const double & incoming_)
  : MultiRegionExactSolution<dim>(
      interface_positions_, source_, sigma_, direction_, incoming_)
{
}

/**
 * \brief Adjust distances travelled in region for transient problem as the
 *        reference solution value may have been inside the domain, in which
 *        case the distances between the boundary and reference point
 *        \f$\mathbf{x} - c\mathbf{\Omega}\f$ need to be subtracted.
 *
 * \param[in] s  distances travelled in each region
 */
template <int dim>
void TransientMultiRegionExactSolution<dim>::adjust_distances(
  std::vector<double> & s) const
{
  // get time
  double t = this->get_time();

  // compute sum of distances
  double s_sum = 0.0;
  for (unsigned int i = 0; i < this->Nr; ++i)
    s_sum += s[i];
  // compute excess distance
  double excess_distance = s_sum - this->c * t;
  // if there is an excess, adjust distances
  if (excess_distance > 0.0)
  {
    for (unsigned int i = 0; i < this->Nr; ++i)
    {
      if (excess_distance > s[i])
      {
        excess_distance -= s[i];
        s[i] = 0.0;
      }
      else
      {
        s[i] -= excess_distance;
        excess_distance = 0.0;
        break;
      }
    }
  }
}

/**
 * \brief Determines the reference solution, which will either be the boundary
 *        value or zero.
 *
 * \param[in] p  point at which exact solution is being evaluated
 *
 * \return reference solution, either the boundary value or zero
 */
template <int dim>
double TransientMultiRegionExactSolution<dim>::compute_reference_solution(
  const Point<dim> & p) const
{
  // get time
  double t = this->get_time();

  // determine if x-c*Omega*t lies within domain; if not,
  // then solution hasn't reached x yet
  bool not_in_domain = false;
  for (unsigned int d = 0; d < dim; ++d)
    if (p[d] - this->direction[d] * this->c * t < 0.0)
      not_in_domain = true;

  // compute solution component due to boundary flux
  // initialize to reference solution
  double ub;
  if (not_in_domain)
    ub = this->incoming;
  else
    ub = 0.0;

  return ub;
}
