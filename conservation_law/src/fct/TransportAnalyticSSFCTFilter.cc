/**
 * \file TransportAnalyticSSFCTFilter.cc
 * \brief Provides the function definitions for the
 *        TransportAnalyticSSFCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] problem_parameters_  problem parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] cell_quadrature_  cell quadrature
 * \param[in] limiter_  limiter
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 */
template <int dim>
TransportAnalyticSSFCTFilter<dim>::TransportAnalyticSSFCTFilter(
  const RunParameters & run_parameters_,
  TransportProblemParameters<dim> & problem_parameters_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const QGauss<dim> & cell_quadrature_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const std::map<unsigned int, double> & dirichlet_values_)
  : SteadyStateFCTFilter<dim>(
      run_parameters_, limiter_, dof_handler_, fe_, dirichlet_values_),
    analytic_bounds(problem_parameters_, dof_handler_, fe_, cell_quadrature_)
{
}

/**
 * \brief Filters antidiffusive fluxes.
 *
 * \param[in] solution  solution \f$\mathbf{U}\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}\f$
 * \param[in] cumulative_antidiffusion  cumulative antidiffusion vector
 *            \f$\bar{\mathbf{p}}\f$
 * \param[inout] limiter_matrix  limiter matrix \f$\mathbf{L}\f$
 * \param[inout] antidiffusion_matrix  antidiffusion matrix \f$\mathbf{P}\f$
 */
template <int dim>
void TransportAnalyticSSFCTFilter<dim>::filter_antidiffusive_fluxes(
  const Vector<double> & solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs,
  const Vector<double> & cumulative_antidiffusion,
  SparseMatrix<double> & limiter_matrix,
  SparseMatrix<double> & antidiffusion_matrix)
{
  // compute solution bounds W- and W+
  compute_solution_bounds(solution, low_order_ss_matrix, ss_rhs);

  // compute antidiffusion bounds Q- and Q+
  this->compute_antidiffusion_bounds(analytic_bounds,
                                     solution,
                                     low_order_ss_matrix,
                                     ss_rhs,
                                     cumulative_antidiffusion);

  // limit antidiffusion fluxes
  this->limiter->compute_limiter_matrix(
    antidiffusion_matrix, this->antidiffusion_bounds, limiter_matrix);
  this->limiter->apply_limiter_matrix(limiter_matrix, antidiffusion_matrix);
}

/**
 * \brief Checks to see if the FCT bounds were satisfied.
 *
 * \param[in] new_solution  new solution vector \f$\mathbf{U}^{n+1}\f$
 *
 * \return flag that FCT bounds were satisfied for all filters
 */
template <int dim>
bool TransportAnalyticSSFCTFilter<dim>::check_bounds(
  const Vector<double> & new_solution)
{
  return analytic_bounds.check_bounds(new_solution);
}

/**
 * \brief Gets the lower solution bound.
 *
 * \return lower solution bound vector
 */
template <int dim>
Vector<double> TransportAnalyticSSFCTFilter<dim>::get_lower_solution_bound() const
{
  return analytic_bounds.lower;
}

/**
 * \brief Gets the upper solution bound.
 *
 * \return upper solution bound vector
 */
template <int dim>
Vector<double> TransportAnalyticSSFCTFilter<dim>::get_upper_solution_bound() const
{
  return analytic_bounds.upper;
}

/**
 * \brief Computes bounds to be imposed on the FCT solution, \f$W_i^-\f$ and
 *        \f$W_i^+\f$.
 *
 * \param[in] solution  solution vector \f$\mathbf{U}\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 */
template <int dim>
void TransportAnalyticSSFCTFilter<dim>::compute_solution_bounds(
  const Vector<double> & solution,
  const SparseMatrix<double> &,
  const Vector<double> &)
{
  // compute analytic solution bounds
  analytic_bounds.update(solution, 0.0, 0.0);
}
