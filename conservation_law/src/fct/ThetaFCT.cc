/**
 * \file ThetaFCT.cc
 * \brief Provides the function definitions for the ThetaFCT class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] consistent_mass_matrix_  consistent mass matrix
 *   \f$\mathbf{M}^C\f$
 * \param[in] lumped_mass_matrix_  lumped mass matrix
 *   \f$\mathbf{M}^L\f$
 */
template <int dim>
ThetaFCT<dim>::ThetaFCT(const RunParameters & run_parameters_,
                        const DoFHandler<dim> & dof_handler_,
                        const FESystem<dim> & fe_,
                        const SparseMatrix<double> & consistent_mass_matrix_,
                        const SparseMatrix<double> & lumped_mass_matrix_)
  : FCT<dim>(run_parameters_, dof_handler_, fe_),
    consistent_mass_matrix(&consistent_mass_matrix_),
    lumped_mass_matrix(&lumped_mass_matrix_),
    theta(run_parameters_.theta),
    use_cumulative_antidiffusion_algorithm(
      run_parameters_.use_cumulative_antidiffusion_algorithm)
{
  // resize temporary vector
  tmp_vector.reinit(this->n_dofs);
}

/**
 * \brief Computes the antidiffusion matrix \f$\mathbf{P}\f$.
 *
 * \param[in] high_order_solution  high-order solution \f$\mathbf{U}^H\f$
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] low_order_diffusion_matrix  new low-order diffusion matrix
 *   \f$\mathbf{D}^{L,n+1}\f$
 * \param[in] low_order_diffusion_matrix  old low-order diffusion matrix
 *   \f$\mathbf{D}^{L,n}\f$
 * \param[in] high_order_diffusion_matrix_new  new high-order diffusion matrix
 *   \f$\mathbf{D}^{H,n+1}\f$
 * \param[in] high_order_diffusion_matrix_old  old high-order diffusion matrix
 *   \f$\mathbf{D}^{H,n}\f$
 * \param[out] antidiffusion_matrix  antidiffusion matrix \f$\mathbf{P}\f$
 */
template <int dim>
void ThetaFCT<dim>::compute_antidiffusion_matrix(
  const Vector<double> & high_order_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const SparseMatrix<double> & low_order_diffusion_matrix_new,
  const SparseMatrix<double> & low_order_diffusion_matrix_old,
  const SparseMatrix<double> & high_order_diffusion_matrix_new,
  const SparseMatrix<double> & high_order_diffusion_matrix_old)
{
  // reset antidiffusion matrix to zero
  this->antidiffusion_matrix = 0;

  // compute time derivative of high-order solution
  Vector<double> & dUdt = tmp_vector;
  dUdt = 0;
  dUdt.add(1.0 / dt, high_order_solution, -1.0 / dt, old_solution);

  // iterate over sparse matrix entries
  SparseMatrix<double>::const_iterator it_mass = consistent_mass_matrix->begin();
  SparseMatrix<double>::const_iterator it_low_new =
    low_order_diffusion_matrix_new.begin();
  SparseMatrix<double>::const_iterator it_low_old =
    low_order_diffusion_matrix_old.begin();
  SparseMatrix<double>::const_iterator it_high_new =
    high_order_diffusion_matrix_new.begin();
  SparseMatrix<double>::const_iterator it_high_old =
    high_order_diffusion_matrix_old.begin();
  SparseMatrix<double>::iterator it_flux = this->antidiffusion_matrix.begin();
  SparseMatrix<double>::iterator it_end = this->antidiffusion_matrix.end();

  for (; it_flux != it_end; ++it_flux,
                            ++it_mass,
                            ++it_low_new,
                            ++it_low_old,
                            ++it_high_new,
                            ++it_high_old)
  {
    // get row and column indices
    const unsigned int i = it_flux->row();
    const unsigned int j = it_flux->column();

    // get values
    const double Mij = it_mass->value();
    const double DLij_new = it_low_new->value();
    const double DLij_old = it_low_old->value();
    const double DHij_new = it_high_new->value();
    const double DHij_old = it_high_old->value();

    // compute flux correction entry
    const double Pij = -Mij * (dUdt(j) - dUdt(i)) +
      (1.0 - theta) * (DLij_old - DHij_old) *
        (old_solution(j) - old_solution(i)) +
      theta * (DLij_new - DHij_new) *
        (high_order_solution(j) - high_order_solution(i));

    // store value
    it_flux->value() = Pij;
  }
}

/**
 * \brief Computes the limited antidiffusion vector with entries
 * \f$\bar{p}_i=\sum\limits_j L_{i,j}P_{i,j}\f$.
 *
 * \param[in] new_solution  new solution \f$\mathbf{U}^{n+1}\f$
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs_new  new steady-state right hand side vector
 *            \f$\mathbf{b}^{n+1}\f$
 * \param[in] ss_rhs_old  old steady-state right hand side vector
 *            \f$\mathbf{b}^n\f$
 * \param[out] antidiffusion_vector  the antidiffusion vector
 *             \f$\bar{\mathbf{p}}\f$
 */
template <int dim>
void ThetaFCT<dim>::compute_antidiffusion_vector(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs_new,
  const Vector<double> & ss_rhs_old,
  Vector<double> & antidiffusion_vector)
{
  // initialize limited antidiffusion matrix
  this->limited_antidiffusion_matrix.copy_from(this->antidiffusion_matrix);

  // filter the antidiffusive fluxes
  filter_antidiffusive_fluxes(new_solution,
                              old_solution,
                              dt,
                              low_order_ss_matrix,
                              ss_rhs_new,
                              ss_rhs_old,
                              this->cumulative_antidiffusion,
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
 * \brief Checks to see if the FCT bounds were satisfied.
 *
 * \param[in] new_solution  new solution vector \f$\mathbf{U}^{n+1}\f$
 *
 * \return flag that FCT bounds were satisfied for all filters
 */
template <int dim>
bool ThetaFCT<dim>::check_bounds(const Vector<double> & new_solution) const
{
  // loop over filters
  bool bounds_satisfied_all_filters = true;
  for (unsigned int k = 0; k < this->n_filters; ++k)
  {
    const bool bounds_satisfied = filters[k]->check_bounds(new_solution);
    bounds_satisfied_all_filters =
      bounds_satisfied_all_filters && bounds_satisfied;
  }

  // return boolean for satisfaction of FCT bounds
  return bounds_satisfied_all_filters;
}

/**
 * \brief Returns lower solution bound vector.
 *
 * \return lower solution bound vector
 */
template <int dim>
Vector<double> ThetaFCT<dim>::get_lower_solution_bound() const
{
  // make sure there is only one filter so there is no amibuity of which
  // bounds to output
  Assert(this->n_filters == 1, ExcNotImplemented());

  return filters[0]->get_lower_solution_bound();
}

/**
 * \brief Returns upper solution bound vector.
 *
 * \return upper solution bound vector
 */
template <int dim>
Vector<double> ThetaFCT<dim>::get_upper_solution_bound() const
{
  // make sure there is only one filter so there is no amibuity of which
  // bounds to output
  Assert(this->n_filters == 1, ExcNotImplemented());

  return filters[0]->get_upper_solution_bound();
}

/**
 * \brief Creates a vector of FCT filters.
 */
template <int dim>
void ThetaFCT<dim>::create_filters()
{
  // create a vector of filter identifier strings
  std::vector<std::string> filter_strings =
    this->create_filter_string_list(this->filter_sequence_string);

  // create vector of filters
  for (unsigned int k = 0; k < this->n_filters; ++k)
  {
    std::shared_ptr<ThetaFCTFilter<dim>> filter =
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
std::shared_ptr<ThetaFCTFilter<dim>> ThetaFCT<dim>::create_filter(
  const std::string & filter_string)
{
  // create filter
  std::shared_ptr<ThetaFCTFilter<dim>> filter;
  if (filter_string == "dmp")
    filter = std::make_shared<DMPThetaFCTFilter<dim>>(
      this->limiter, *this->dof_handler, *this->fe, *lumped_mass_matrix, theta);
  else
    throw ExcNotImplemented();

  // return pointer to new filter
  return filter;
}

/**
 * \brief Filters the antidiffusive fluxes \f$P_{i,j}\f$.
 *
 * \param[in] new_solution  new solution \f$\mathbf{U}^{n+1}\f$
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs_new  new steady-state right hand side vector
 *            \f$\mathbf{b}^{n+1}\f$
 * \param[in] ss_rhs_old  old steady-state right hand side vector
 *            \f$\mathbf{b}^n\f$
 * \param[in] cumulative_antidiffusion  cumulative antidiffusion vector
 *            \f$\bar{\mathbf{p}}\f$
 * \param[inout] limiter_matrix  limiter matrix \f$\mathbf{L}\f$
 * \param[inout] antidiffusion_matrix  antidiffusion matrix \f$\mathbf{P}\f$
 */
template <int dim>
void ThetaFCT<dim>::filter_antidiffusive_fluxes(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs_new,
  const Vector<double> & ss_rhs_old,
  const Vector<double> & cumulative_antidiffusion,
  SparseMatrix<double> & limiter_matrix,
  SparseMatrix<double> & antidiffusion_matrix)
{
  // loop over filters
  for (unsigned int k = 0; k < this->n_filters; ++k)
    filters[k]->filter_antidiffusive_fluxes(new_solution,
                                            old_solution,
                                            dt,
                                            low_order_ss_matrix,
                                            ss_rhs_new,
                                            ss_rhs_old,
                                            cumulative_antidiffusion,
                                            limiter_matrix,
                                            antidiffusion_matrix);
}
