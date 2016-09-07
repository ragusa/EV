/**
 * \file TransportSteadyStateFCT.cc
 * \brief Provides the function definitions for the TransportSteadyStateFCT
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
 * \param[in] cell_quadrature_  cell quadrature
 * \param[in] dx_min_  minimum cell diameter
 */
template <int dim>
TransportSteadyStateFCT<dim>::TransportSteadyStateFCT(
  const TransportRunParameters & run_parameters_,
  TransportProblemParameters<dim> & problem_parameters_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const std::map<unsigned int, double> & dirichlet_values_,
  const QGauss<dim> & cell_quadrature_,
  const double & dx_min_)
  : SteadyStateFCT<dim>(run_parameters_, dof_handler_, fe_, dirichlet_values_),
    problem_parameters(&problem_parameters_),
    cell_quadrature(&cell_quadrature_),
    dx_min(dx_min_)
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
std::shared_ptr<SteadyStateFCTFilter<dim>> TransportSteadyStateFCT<
  dim>::create_filter(const std::string & filter_string)
{
  // create filter
  std::shared_ptr<SteadyStateFCTFilter<dim>> filter;
  if (filter_string == "dmp")
    filter =
      std::make_shared<DMPSteadyStateFCTFilter<dim>>(*this->run_parameters,
                                                     this->limiter,
                                                     *this->dof_handler,
                                                     *this->fe,
                                                     *this->dirichlet_values);
  else if (filter_string == "dmp_analytic")
    filter = std::make_shared<TransportDMPAnalyticSSFCTFilter<dim>>(
      *this->run_parameters,
      *problem_parameters,
      *this->dof_handler,
      *this->fe,
      *cell_quadrature,
      this->limiter,
      *this->dirichlet_values);
  else if (filter_string == "analytic")
    filter = std::make_shared<TransportAnalyticSSFCTFilter<dim>>(
      *this->run_parameters,
      *problem_parameters,
      *this->dof_handler,
      *this->fe,
      *cell_quadrature,
      this->limiter,
      *this->dirichlet_values,
      dx_min);
  else if (filter_string == "upwind_analytic")
    filter = std::make_shared<TransportUpwindAnalyticSSFCTFilter<dim>>(
      *this->run_parameters,
      *problem_parameters,
      *this->dof_handler,
      *this->fe,
      this->limiter,
      *this->dirichlet_values,
      dx_min);
  else if (filter_string == "dmp_exact")
    filter =
      std::make_shared<TransportExactDMPSSFCTFilter<dim>>(*this->run_parameters,
                                                          this->limiter,
                                                          *this->dof_handler,
                                                          *this->fe,
                                                          *this->dirichlet_values,
                                                          *problem_parameters);
  else if (filter_string == "analytic_exact")
    filter = std::make_shared<TransportExactAnalyticSSFCTFilter<dim>>(
      *this->run_parameters,
      *problem_parameters,
      *this->dof_handler,
      *this->fe,
      *cell_quadrature,
      this->limiter,
      *this->dirichlet_values,
      dx_min);
  else
    AssertThrow(false, ExcNotImplemented());

  // return pointer to new filter
  return filter;
}
