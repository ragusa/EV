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

/**
 * \brief Returns the size of the viscosity map
 *
 * \return size of the viscosity map
 */
template <int dim>
unsigned int Viscosity<dim>::size() const
{
  return values.size();
}

/**
 * \brief Returns the cell map of viscosity values
 *
 * \return cell map of viscosity values
 */
template <int dim>
typename Viscosity<dim>::CellMap Viscosity<dim>::get_values() const
{
  return values;
}
