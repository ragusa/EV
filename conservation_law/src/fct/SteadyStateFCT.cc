/**
 * \file SteadyStateFCT.cc
 * \brief Provides the function definitions for the SteadyStateFCT class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 */
template <int dim>
SteadyStateFCT<dim>::SteadyStateFCT(const RunParameters & run_parameters_,
                                    const DoFHandler<dim> & dof_handler_,
                                    const FESystem<dim> & fe_)
  : FCT<dim>(run_parameters_, dof_handler_, fe_),
    use_cumulative_antidiffusion_algorithm(
      run_parameters_.use_cumulative_antidiffusion_algorithm)
{
}

/**
 * \brief Computes the antidiffusion matrix \f$\mathbf{P}\f$.
 *
 * \param[in] high_order_solution  high-order solution \f$\mathbf{U}^H\f$
 * \param[in] low_order_diffusion_matrix  low-order diffusion matrix
 *   \f$\mathbf{D}^L\f$
 * \param[in] high_order_diffusion_matrix  high-order diffusion matrix
 *   \f$\mathbf{D}^H\f$
 */
template <int dim>
void SteadyStateFCT<dim>::compute_antidiffusion_matrix(
  const Vector<double> & high_order_solution,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const SparseMatrix<double> & high_order_diffusion_matrix)
{
  // reset antidiffusion matrix to zero
  this->antidiffusion_matrix = 0;

  // iterate over sparse matrix entries
  SparseMatrix<double>::const_iterator it_low =
    low_order_diffusion_matrix.begin();
  SparseMatrix<double>::const_iterator it_high =
    high_order_diffusion_matrix.begin();
  SparseMatrix<double>::iterator it_flux = this->antidiffusion_matrix.begin();
  SparseMatrix<double>::iterator it_end = this->antidiffusion_matrix.end();

  for (; it_flux != it_end; ++it_flux, ++it_low, ++it_high)
  {
    // get row and column indices
    const unsigned int i = it_flux->row();
    const unsigned int j = it_flux->column();

    // get values
    const double DLij = it_low->value();
    const double DHij = it_high->value();

    // compute flux correction entry
    const double Fij =
      (DLij - DHij) * (high_order_solution(j) - high_order_solution(i));

    // store value
    it_flux->value() = Fij;
  }

  // also reset cumulative antidiffusion vector
  this->cumulative_antidiffusion = 0;
}

/**
 * \brief Computes the limited antidiffusion vector with entries
 * \f$\bar{p}_i=\sum\limits_j L_{i,j}P_{i,j}\f$.
 *
 * \param[in] solution  solution \f$\mathbf{U}\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}\f$
 * \param[out] antidiffusion_vector  the antidiffusion vector
 *             \f$\bar{\mathbf{p}}\f$
 */
template <int dim>
void SteadyStateFCT<dim>::compute_antidiffusion_vector(
  const Vector<double> & solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs,
  Vector<double> & antidiffusion_vector)
{
  // initialize limited antidiffusion matrix
  this->limited_antidiffusion_matrix.copy_from(this->antidiffusion_matrix);

  // filter the antidiffusive fluxes
  filter_antidiffusive_fluxes(solution,
                              low_order_ss_matrix,
                              ss_rhs,
                              this->limiter_matrix,
                              this->limited_antidiffusion_matrix);

  // compute sum of limited antidiffusive fluxes
  this->compute_row_sum_vector(this->limited_antidiffusion_matrix,
                               antidiffusion_vector);

  // add cumulative antidiffusion of previous iterations if requested
  if (use_cumulative_antidiffusion_algorithm)
  {
    antidiffusion_vector.add(1.0, this->cumulative_antidiffusion);
    this->cumulative_antidiffusion = antidiffusion_vector;
    this->antidiffusion_matrix.add(-1.0, this->limited_antidiffusion_matrix);
  }
}

/**
 * \brief Creates a vector of FCT filters.
 */
template <int dim>
void SteadyStateFCT<dim>::create_filters()
{
  // create a vector of filter identifier strings
  std::vector<std::string> filter_strings =
    this->create_filter_string_list(this->filter_sequence_string);

  // create vector of filters
  for (unsigned int k = 0; k < this->n_filters; ++k)
  {
    std::shared_ptr<SteadyStateFCTFilter<dim>> filter =
      create_filter(filter_strings[k]);
    filters.push_back(filter);
  }
}

/**
 * \brief Creates an FCT filter.
 *
 * \param[in] filter_string  string identifier for filter to create
 *
 * \return pointer to created FCT filter
 */
template <int dim>
std::shared_ptr<SteadyStateFCTFilter<dim>> SteadyStateFCT<dim>::create_filter(
  const std::string & filter_string)
{
  // create filter
  std::shared_ptr<SteadyStateFCTFilter<dim>> filter;
  if (filter_string == "dmp")
    filter = std::make_shared<DMPSteadyStateFCTFilter<dim>>(
      this->limiter, *this->dof_handler, *this->fe);
  else
    throw ExcNotImplemented();

  // return pointer to new filter
  return filter;
}

/**
 * \brief Filters the antidiffusive fluxes \f$P_{i,j}\f$.
 *
 * \param[in] solution  solution \f$\mathbf{U}\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}\f$
 * \param[inout] limiter_matrix  limiter matrix \f$\mathbf{L}\f$
 * \param[inout] antidiffusion_matrix  antidiffusion matrix \f$\mathbf{P}\f$
 */
template <int dim>
void SteadyStateFCT<dim>::filter_antidiffusive_fluxes(
  const Vector<double> & solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs,
  SparseMatrix<double> & limiter_matrix,
  SparseMatrix<double> & antidiffusion_matrix)
{
  // loop over filters
  for (unsigned int k = 0; k < this->n_filters; ++k)
    filters[k]->filter_antidiffusive_fluxes(solution,
                                            low_order_ss_matrix,
                                            ss_rhs,
                                            this->cumulative_antidiffusion,
                                            limiter_matrix,
                                            antidiffusion_matrix);
}
