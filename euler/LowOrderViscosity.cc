/**
 * \file LowOrderViscosity.cc
 * \brief Provides the function definitions for the LowOrderViscosity class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] c_max_ coefficient for low-order viscosity
 * \param[in] cell_diameter_ cell diameter for each cell
 * \param[in] max_flux_speed_ max flux speed for each cell
 * \param[in] dof_handler_ degree of freedom handler
 */
template <int dim>
LowOrderViscosity<dim>::LowOrderViscosity(const double & c_max_,
                                          CellMap & cell_diameter_,
                                          CellMap & max_flux_speed_,
                                          const DoFHandler<dim> & dof_handler_)
  : Viscosity<dim>(dof_handler_),
    c_max(c_max_),
    cell_diameter(&cell_diameter_),
    max_flux_speed(&max_flux_speed_)
{
}

template <int dim>
void LowOrderViscosity<dim>::update()
{
  // loop over cells to compute low-order viscosity at each quadrature point
  Cell cell = this->dof_handler->begin_active(), endc = this->dof_handler->end();
  for (; cell != endc; ++cell)
  {
    this->values[cell] =
      std::abs(c_max * (*cell_diameter)[cell] * (*max_flux_speed)[cell]);
  }
}
