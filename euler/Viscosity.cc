/**
 * \file Viscosity.cc
 * \brief Provides the function definitions for the Viscosity class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] dof_handler_ degree of freedom handler
 */
template <int dim>
Viscosity<dim>::Viscosity(const DoFHandler<dim> & dof_handler_)
  : dof_handler(&dof_handler_)
{
}

/**
 * \brief Accesses the viscosity value for the given cell
 *
 * \param[in] cell cell iterator
 *
 * \return viscosity value for the cell
 */
template <int dim>
double & Viscosity<dim>::operator[](const Cell & cell)
{
  return values[cell];
}
