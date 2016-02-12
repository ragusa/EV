/**
 * \file ShallowWaterViscosityMultiplier.cc
 * \brief Provides the function definitions for the
 * ShallowWaterViscosityMultiplier class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ShallowWaterViscosityMultiplier<dim>::ShallowWaterViscosityMultiplier(
  const double & gravity_)
  : ViscosityMultiplier<dim>(),
    gravity(gravity_),
    height_extractor(0),
    momentum_extractor(1)
{
}

/**
 * \brief Computes max Froude number in cell.
 *
 * The max Froude number is taken to be the max over all quadrature points
 * in the cell:
 * \f[
 *   \mbox{Fr}_K = \max\limits_q \frac{u_q}{a_q} \,.
 * \f]
 */
template <int dim>
double ShallowWaterViscosityMultiplier<dim>::get_multiplier(
  const FEValues<dim> & fe_values, const Vector<double> & solution) const
{
  // get number of quadrature points
  const unsigned int n_quadrature_points = fe_values.n_quadrature_points;

  // extract height
  std::vector<double> height(n_quadrature_points);
  fe_values[height_extractor].get_function_values(solution, height);

  // extract momentum
  std::vector<Tensor<1, dim>> momentum(n_quadrature_points);
  fe_values[momentum_extractor].get_function_values(solution, momentum);

  // compute maximum Froude number for cell
  double froude_cell = 0.0;
  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    // compute Froude number at quadrature point
    const Tensor<1, dim> velocity = momentum[q] / height[q];
    const double speed = velocity.norm();
    const double speed_of_sound = 0.5 * gravity * height[q];
    const double froude = speed / speed_of_sound;

    // take maximum over all quadrature points
    froude_cell = std::max(froude_cell, froude);
  }

  return froude_cell;
}
