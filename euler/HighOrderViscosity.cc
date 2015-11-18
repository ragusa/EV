/**
 * \file HighOrderViscosity.cc
 * \brief Provides the function definitions for the HighOrderViscosity class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] low_order_viscosity_ pointer to low-order viscosity
 * \param[in] entropy_viscosity_ pointer to entropy viscosity
 * \param[in] use_low_order_for_first_step_ flag to use low-order viscosity
 *            for first time step
 * \param[in] dof_handler_ degree of freedom handler
 */
template <int dim>
HighOrderViscosity<dim>::HighOrderViscosity(
  const std::shared_ptr<Viscosity<dim>> & low_order_viscosity_,
  const std::shared_ptr<Viscosity<dim>> & entropy_viscosity_,
  const bool & use_low_order_for_first_step_,
  const DoFHandler<dim> & dof_handler_)
  : Viscosity<dim>(dof_handler_),
    low_order_viscosity(low_order_viscosity_),
    entropy_viscosity(entropy_viscosity_),
    use_low_order_for_first_step(use_low_order_for_first_step_)
{
}

/**
 * \brief Computes high-order viscosity values.
 *
 * \param[in] new_solution new solution
 * \param[in] old_solution old solution
 * \param[in] dt time step size
 * \param[in] n time index
 */
template <int dim>
void HighOrderViscosity<dim>::update(const Vector<double> & new_solution,
                                   const Vector<double> & old_solution,
                                   const double & dt,
const unsigned int & n,
{
  if (use_low_order_for_first_step && n == 1)
  {
    // update the low-order viscosity
    low_order_viscosity->update(new_solution, old_solution, dt, n);

    // copy the low-order viscosity values
    Cell cell = this->dof_handler->begin_active(),
         endc = this->dof_handler->end();
    for (; cell != endc; ++cell)
      this->values[cell] = (*low_order_viscosity)[cell];
  }
  else
  {
    // update the low-order viscosity
    low_order_viscosity->update(new_solution, old_solution, dt, n);

    // update the entropy viscosity
    entropy_viscosity->update(new_solution, old_solution, dt, n);

    // take the minimum of the low-order viscosity and entropy viscosity
    Cell cell = this->dof_handler->begin_active(),
         endc = this->dof_handler->end();
    for (; cell != endc; ++cell)
      this->values[cell] =
        std::min((*low_order_viscosity)[cell], (*entropy_viscosity)[cell]);
  }
}
