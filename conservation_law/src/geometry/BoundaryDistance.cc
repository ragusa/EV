/**
 * \file BoundaryDistance.cc
 * \brief Provides the function definitions for the BoundaryDistance class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] direction_  transport direction vector
 */
template <int dim>
BoundaryDistance<dim>::BoundaryDistance(const Tensor<1, dim> & direction_)
  : direction(direction_)
{
  // assert that all components of direction are non-negative
  for (unsigned int d = 0; d < dim; ++d)
  {
    AssertThrow(direction[d] >= 0.0, ExcAssumptionViolated());
  }
}
