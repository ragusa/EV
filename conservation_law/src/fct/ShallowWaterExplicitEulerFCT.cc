/**
 * \file ShallowWaterExplicitEulerFCT.cc
 * \brief Provides the function definitions for the ShallowWaterExplicitEulerFCT
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] problem_parameters_  problem parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 * \param[in] consistent_mass_matrix_  consistent mass matrix
 *   \f$\mathbf{M}^C\f$
 * \param[in] lumped_mass_matrix_  lumped mass matrix
 *   \f$\mathbf{M}^L\f$
 */
template <int dim>
ShallowWaterExplicitEulerFCT<dim>::ShallowWaterExplicitEulerFCT(
  const ShallowWaterRunParameters & run_parameters_,
  const ShallowWaterProblemParameters<dim> & problem_parameters_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const std::map<unsigned int, double> & dirichlet_values_,
  const SparseMatrix<double> & consistent_mass_matrix_,
  const SparseMatrix<double> & lumped_mass_matrix_)
  : ExplicitEulerFCT<dim>(run_parameters_,
                          dof_handler_,
                          fe_,
                          dirichlet_values_,
                          consistent_mass_matrix_,
                          lumped_mass_matrix_),
    gravity(problem_parameters_.gravity)
{
  // create FCT filters
  this->create_filters();
}

/**
 * \brief Creates an FCT filter.
 *
 * \param[in] filter_string  string identifier for filter to create
 *
 * \return pointer to created FCT filter
 */
template <int dim>
std::shared_ptr<ExplicitEulerFCTFilter<dim>> ShallowWaterExplicitEulerFCT<
  dim>::create_filter(const std::string & filter_string)
{
  // create filter
  std::shared_ptr<ExplicitEulerFCTFilter<dim>> filter;
  if (filter_string == "characteristic")
    filter =
      std::make_shared<SWCharacteristicFCTFilter<dim>>(*this->run_parameters,
                                                       this->limiter,
                                                       *this->dof_handler,
                                                       *this->fe,
                                                       *this->lumped_mass_matrix,
                                                       gravity,
                                                       *this->dirichlet_values);
  else
    AssertThrow(false, ExcNotImplemented());

  // return pointer to new filter
  return filter;
}
