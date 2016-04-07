/**
 * \file TransportExplicitEulerFCT.cc
 * \brief Provides the function definitions for the TransportExplicitEulerFCT
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] problem_parameters_  problem parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] consistent_mass_matrix_  consistent mass matrix
 *   \f$\mathbf{M}^C\f$
 * \param[in] lumped_mass_matrix_  lumped mass matrix
 *   \f$\mathbf{M}^L\f$
 */
template <int dim>
TransportExplicitEulerFCT<dim>::TransportExplicitEulerFCT(
  const TransportRunParameters & run_parameters_,
  const TransportProblemParameters<dim> &,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const SparseMatrix<double> & consistent_mass_matrix_,
  const SparseMatrix<double> & lumped_mass_matrix_)
  : ExplicitEulerFCT<dim>(run_parameters_,
                          dof_handler_,
                          fe_,
                          consistent_mass_matrix_,
                          lumped_mass_matrix_)
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
std::shared_ptr<ExplicitEulerFCTFilter<dim>> TransportExplicitEulerFCT<
  dim>::create_filter(const std::string & filter_string)
{
  // create filter
  std::shared_ptr<ExplicitEulerFCTFilter<dim>> filter;
  if (filter_string == "dmp")
    filter =
      std::make_shared<DMPExplicitEulerFCTFilter<dim>>(*this->run_parameters,
                                                       this->limiter,
                                                       *this->dof_handler,
                                                       *this->fe,
                                                       *this->lumped_mass_matrix);
  else
    AssertThrow(false, ExcNotImplemented());

  // return pointer to new filter
  return filter;
}
