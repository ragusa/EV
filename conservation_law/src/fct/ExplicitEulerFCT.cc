/**
 * \file ExplicitEulerFCT.cc
 * \brief Provides the function definitions for the ExplicitEulerFCT class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 * \param[in] consistent_mass_matrix_  consistent mass matrix
 *   \f$\mathbf{M}^C\f$
 * \param[in] lumped_mass_matrix_  lumped mass matrix
 *   \f$\mathbf{M}^L\f$
 */
template <int dim>
ExplicitEulerFCT<dim>::ExplicitEulerFCT(
  const RunParameters & run_parameters_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const std::map<unsigned int, double> & dirichlet_values_,
  const SparseMatrix<double> & consistent_mass_matrix_,
  const SparseMatrix<double> & lumped_mass_matrix_)
  : FCT<dim>(run_parameters_, dof_handler_, fe_, dirichlet_values_),
    consistent_mass_matrix(&consistent_mass_matrix_),
    lumped_mass_matrix(&lumped_mass_matrix_)
{
  // resize temporary vector
  tmp_vector.reinit(this->n_dofs);
}

/**
 * \brief Computes the limited antidiffusion vector with entries
 * \f$\bar{p}_i=\sum\limits_j L_{i,j}P_{i,j}\f$.
 *
 * \param[in] high_order_solution  high-order solution \f$\mathbf{U}^H\f$
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] inviscid_ss_flux  inviscid steady-state flux vector
 *            (entries are \f$(\mathbf{A}\mathbf{U}^n)_i\f$ for scalar case,
 *            \f$\sum\limits_j\mathbf{c}_{i,j}\cdot\mathrm{F}^n_j\f$ for
 *            systems case)
 * \param[in] ss_reaction  steady-state reaction vector \f$\mathbf{\sigma}\f$,
 *            where \f$\sigma_i = \int\limits_{S_i}
 *            \varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *            = \sum\limits_j\int\limits_{S_{i,j}}
 *            \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \f$
 * \param[in] low_order_diffusion_matrix  low-order diffusion matrix
 *            \f$\mathbf{D}^L\f$
 * \param[in] high_order_diffusion_matrix  high-order diffusion matrix
 *            \f$\mathbf{D}^H\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}^n\f$
 * \param[out] antidiffusion_vector  the antidiffusion vector
 *             \f$\bar{\mathbf{p}}\f$
 */
template <int dim>
void ExplicitEulerFCT<dim>::compute_antidiffusion_vector(
  const Vector<double> & high_order_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> & inviscid_ss_flux,
  const Vector<double> & ss_reaction,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const SparseMatrix<double> & high_order_diffusion_matrix,
  const Vector<double> & ss_rhs,
  Vector<double> & antidiffusion_vector)
{
  // compute antidiffusive fluxes
  compute_antidiffusion_matrix(high_order_solution,
                               old_solution,
                               dt,
                               low_order_diffusion_matrix,
                               high_order_diffusion_matrix,
                               this->antidiffusion_matrix);

  // filter the antidiffusive fluxes
  filter_antidiffusive_fluxes(old_solution,
                              dt,
                              inviscid_ss_flux,
                              ss_reaction,
                              low_order_diffusion_matrix,
                              ss_rhs,
                              this->limiter_matrix,
                              this->antidiffusion_matrix);

  // compute sum of limited antidiffusive fluxes
  this->compute_row_sum_vector(this->antidiffusion_matrix, antidiffusion_vector);
}

/**
 * \brief Checks to see if the FCT bounds were satisfied.
 *
 * \param[in] new_solution  new solution vector \f$\mathbf{U}^{n+1}\f$
 *
 * \return flag that FCT bounds were satisfied for all filters
 */
template <int dim>
bool ExplicitEulerFCT<dim>::check_bounds(
  const Vector<double> & new_solution) const
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
Vector<double> ExplicitEulerFCT<dim>::get_lower_solution_bound() const
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
Vector<double> ExplicitEulerFCT<dim>::get_upper_solution_bound() const
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
void ExplicitEulerFCT<dim>::create_filters()
{
  // create a vector of filter identifier strings
  std::vector<std::string> filter_strings =
    this->create_filter_string_list(this->filter_sequence_string);

  // create vector of filters
  for (unsigned int k = 0; k < this->n_filters; ++k)
  {
    std::shared_ptr<ExplicitEulerFCTFilter<dim>> filter =
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
std::shared_ptr<ExplicitEulerFCTFilter<dim>> ExplicitEulerFCT<dim>::create_filter(
  const std::string & filter_string)
{
  // create filter
  std::shared_ptr<ExplicitEulerFCTFilter<dim>> filter;
  if (filter_string == "dmp")
    filter =
      std::make_shared<DMPExplicitEulerFCTFilter<dim>>(*this->run_parameters,
                                                       this->limiter,
                                                       *this->dof_handler,
                                                       *this->fe,
                                                       *lumped_mass_matrix);
  else
    AssertThrow(false, ExcNotImplemented());

  // return pointer to new filter
  return filter;
}

/**
 * \brief Computes the antidiffusion matrix \f$\mathbf{P}\f$.
 *
 * \param[in] high_order_solution  high-order solution \f$\mathbf{U}^H\f$
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] low_order_diffusion_matrix  low-order diffusion matrix
 *   \f$\mathbf{D}^L\f$
 * \param[in] high_order_diffusion_matrix  high-order diffusion matrix
 *   \f$\mathbf{D}^H\f$
 * \param[out] antidiffusion_matrix  antidiffusion matrix \f$\mathbf{P}\f$
 */
template <int dim>
void ExplicitEulerFCT<dim>::compute_antidiffusion_matrix(
  const Vector<double> & high_order_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const SparseMatrix<double> & high_order_diffusion_matrix,
  SparseMatrix<double> & antidiffusion_matrix)
{
  // reset antidiffusion matrix to zero
  antidiffusion_matrix = 0;

  // compute time derivative of high-order solution
  Vector<double> & dUdt = tmp_vector;
  dUdt = 0;
  dUdt.add(1.0 / dt, high_order_solution, -1.0 / dt, old_solution);

  // iterate over sparse matrix entries
  SparseMatrix<double>::const_iterator it_mass = consistent_mass_matrix->begin();
  SparseMatrix<double>::const_iterator it_low =
    low_order_diffusion_matrix.begin();
  SparseMatrix<double>::const_iterator it_high =
    high_order_diffusion_matrix.begin();
  SparseMatrix<double>::iterator it_flux = antidiffusion_matrix.begin();
  SparseMatrix<double>::iterator it_end = antidiffusion_matrix.end();

  for (; it_flux != it_end; ++it_flux, ++it_mass, ++it_low, ++it_high)
  {
    // get row and column indices
    const unsigned int i = it_flux->row();
    const unsigned int j = it_flux->column();

    // get values
    const double Mij = it_mass->value();
    const double DLij = it_low->value();
    const double DHij = it_high->value();

    // compute flux correction entry
    const double Fij = -Mij * (dUdt(j) - dUdt(i)) +
      (DLij - DHij) * (old_solution(j) - old_solution(i));

    // store value
    it_flux->value() = Fij;
  }
}

/**
 * \brief Filters the antidiffusive fluxes \f$P_{i,j}\f$.
 *
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] inviscid_ss_flux  inviscid steady-state flux vector
 *            (entries are \f$(\mathbf{A}\mathbf{U}^n)_i\f$ for scalar case,
 *            \f$\sum\limits_j\mathbf{c}_{i,j}\cdot\mathrm{F}^n_j\f$ for
 *            systems case)
 * \param[in] ss_reaction  steady-state reaction vector \f$\mathbf{\sigma}\f$,
 *            where \f$\sigma_i = \int\limits_{S_i}
 *            \varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *            = \sum\limits_j\int\limits_{S_{i,j}}
 *            \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \f$
 * \param[in] low_order_diffusion_matrix  low-order diffusion matrix
 *            \f$\mathbf{D}^L\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}^n\f$
 * \param[inout] limiter_matrix  limiter matrix \f$\mathbf{L}\f$
 * \param[inout] antidiffusion_matrix  antidiffusion matrix \f$\mathbf{P}\f$
 */
template <int dim>
void ExplicitEulerFCT<dim>::filter_antidiffusive_fluxes(
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> & inviscid_ss_flux,
  const Vector<double> & ss_reaction,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const Vector<double> & ss_rhs,
  SparseMatrix<double> & limiter_matrix,
  SparseMatrix<double> & antidiffusion_matrix)
{
  // loop over filters
  for (unsigned int k = 0; k < this->n_filters; ++k)
    filters[k]->filter_antidiffusive_fluxes(old_solution,
                                            dt,
                                            inviscid_ss_flux,
                                            ss_reaction,
                                            low_order_diffusion_matrix,
                                            ss_rhs,
                                            limiter_matrix,
                                            antidiffusion_matrix);
}
