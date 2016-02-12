/**
 * \file ViscosityMultiplier.cc
 * \brief Provides the function definitions for the ViscosityMultiplier class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ViscosityMultiplier<dim>::ViscosityMultiplier()
{
}

/**
 * \brief Returns default value of 1 as the multiplier for the cell viscosity.
 *
 * \param[in] fe_values_ FE values
 * \param[in] solution_ solution vector
 */
template <int dim>
double ViscosityMultiplier<dim>::get_multiplier(const FEValues<dim> &,
                                                const Vector<double> &) const
{
  return 1.0;
}
