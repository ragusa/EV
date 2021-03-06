/**
 * \file ConstantViscosity.cc
 * \brief Provides the function definitions for the ConstantViscosity class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] constant_value_ constant value to use for viscosity
 * \param[in] dof_handler_ degree of freedom handler
 */
template <int dim>
ConstantViscosity<dim>::ConstantViscosity(const double & constant_value_,
                                          const DoFHandler<dim> & dof_handler_)
  : Viscosity<dim>(dof_handler_), constant_value(constant_value_)
{
  Cell cell = this->dof_handler->begin_active(), endc = this->dof_handler->end();
  for (; cell != endc; ++cell)
    this->values[cell] = constant_value;
}

/**
 * \brief Does nothing, as values are constant and already set.
 *
 * \param[in] new_solution new solution
 * \param[in] old_solution old solution
 * \param[in] dt time step size
 * \param[in] n time index
 */
template <int dim>
void ConstantViscosity<dim>::update(const Vector<double> &,
                                    const Vector<double> &,
                                    const double &,
                                    const unsigned int &)
{
}
