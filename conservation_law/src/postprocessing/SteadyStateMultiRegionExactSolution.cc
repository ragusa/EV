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
SteadyStateMultiRegionExactSolution<dim>::SteadyStateMultiRegionExactSolution(
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
void SteadyStateMultiRegionExactSolution<dim>::adjust_distances(
  std::vector<double> &) const
{
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
double SteadyStateMultiRegionExactSolution<dim>::compute_reference_solution(
  const Point<dim> &) const
{
  return this->incoming;
}
