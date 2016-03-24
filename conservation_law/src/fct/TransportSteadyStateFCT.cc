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
 * \param[in] cell_quadrature_  cell quadrature
 */
template <int dim>
TransportSteadyStateFCT<dim>::TransportSteadyStateFCT(
  const TransportRunParameters & run_parameters_,
  TransportProblemParameters<dim> & problem_parameters_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const QGauss<dim> & cell_quadrature_)
  : SteadyStateFCT<dim>(run_parameters_, dof_handler_, fe_),
    problem_parameters(&problem_parameters_),
    cell_quadrature(&cell_quadrature_)
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
    filter = std::make_shared<DMPSteadyStateFCTFilter<dim>>(
      *this->run_parameters, this->limiter, *this->dof_handler, *this->fe);
  else if (filter_string == "dmp_analytic")
    filter = std::make_shared<TransportDMPAnalyticSSFCTFilter<dim>>(
      *this->run_parameters,
      *problem_parameters,
      *this->dof_handler,
      *this->fe,
      *cell_quadrature,
      this->limiter);
  else
    throw ExcNotImplemented();

  // return pointer to new filter
  return filter;
}
