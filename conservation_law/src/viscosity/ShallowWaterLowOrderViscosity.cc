/**
 * \file ShallowWaterLowOrderViscosity.cc
 * \brief Provides the function definitions for the
 *        ShallowWaterLowOrderViscosity class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ShallowWaterLowOrderViscosity<dim>::ShallowWaterLowOrderViscosity()
{
}

/**
 * \brief Computes low-order viscosity for each cell.
 *
 * The low-order viscosity is computed as
 * \f[
 *   \nu_K^L = c_{max} h_K \lambda_{K,max} \,,
 * \f]
 * where \f$\lambda_{K,max}\f$ is the maximum flux speed on cell \f$K\f$,
 * or if the user specifies, the low-order viscosity definition above is
 * multiplied by the local Froude number:
 * \f[
 *   \nu_K^L = c_{max} \mbox{Fr}_K h_K \lambda_{K,max} \,,
 * \f]
 * where the Froude number is taken to be the max over all quadrature points:
 * \f[
 *   \mbox{Fr}_K = \max\limits_q \frac{u_q}{a_q} \,.
 * \f]
 *
 * \param[in] using_low_order_scheme flag that the low-order viscosities are
 *            are used in a low-order scheme as opposed to an entropy
 *            viscosity scheme
 */
/*
template <int dim>
void ShallowWater<dim>::update_old_low_order_viscosity(
  const bool & using_low_order_scheme)
{
  const double c_max = this->parameters.first_order_viscosity_coef;

  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();

  if (using_low_order_scheme &&
      sw_parameters.multiply_low_order_viscosity_by_froude)
  {
    FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);

    // compute low-order viscosity for each cell
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      // extract height and momentum
      std::vector<double> height(this->n_q_points_cell);
      fe_values[height_extractor].get_function_values(this->new_solution, height);
      std::vector<Tensor<1, dim>> momentum(this->n_q_points_cell);
      fe_values[momentum_extractor].get_function_values(this->new_solution,
                                                        momentum);

      // compute fluid speed and sound speed
      std::vector<double> speed = compute_speed(height, momentum);
      std::vector<double> sound_speed = compute_sound_speed(height);

      // compute Froude number for cell by taking the maximum over
      // quadrature points
      double froude = 0.0;
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
        froude = std::max(froude, speed[q] / sound_speed[q]);

      this->first_order_viscosity_map[cell] =
        std::abs(c_max * froude * cell->diameter() *
                 this->max_flux_speed_cell[cell]);
    }
  }
  else
  {
    // compute low-order viscosity for each cell
    for (; cell != endc; ++cell)
    {
      this->first_order_viscosity_map[cell] = std::abs(
        c_max * cell->diameter() * this->max_flux_speed_cell[cell]);
    }
  }
}
*/
