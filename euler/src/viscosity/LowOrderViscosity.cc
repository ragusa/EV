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
 * \param[in] fe_ finite element system
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] viscosity_multiplier_ viscosity multiplier
 */
template <int dim>
LowOrderViscosity<dim>::LowOrderViscosity(
  const double & c_max_,
  CellMap & cell_diameter_,
  CellMap & max_flux_speed_,
  const FESystem<dim> & fe_,
  const DoFHandler<dim> & dof_handler_,
  const QGauss<dim> & cell_quadrature_,
  const std::shared_ptr<ViscosityMultiplier<dim>> & viscosity_multiplier_)
  : Viscosity<dim>(dof_handler_),
    c_max(c_max_),
    cell_diameter(&cell_diameter_),
    max_flux_speed(&max_flux_speed_),
    fe(&fe_),
    cell_quadrature(&cell_quadrature_),
    viscosity_multiplier(viscosity_multiplier_)
{
}

template <int dim>
void LowOrderViscosity<dim>::update(const Vector<double> & new_solution,
                                    const Vector<double> &,
                                    const double &,
                                    const unsigned int &)
{
  // FE values
  FEValues<dim> fe_values(*fe, *cell_quadrature, update_values);

  // loop over cells to compute low-order viscosity at each quadrature point
  Cell cell = this->dof_handler->begin_active(), endc = this->dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // get multiplier
    const double multiplier =
      viscosity_multiplier->get_multiplier(fe_values, new_solution);

    // compute viscosity value
    this->values[cell] =
      std::abs(c_max * (*cell_diameter)[cell] * (*max_flux_speed)[cell]) *
      multiplier;
  }
}
