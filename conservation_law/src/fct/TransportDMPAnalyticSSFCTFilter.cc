/**
 * \file TransportDMPAnalyticSSFCTFilter.cc
 * \brief Provides the function definitions for the
 *        TransportDMPAnalyticSSFCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 */
template <int dim>
TransportDMPAnalyticSSFCTFilter<dim>::TransportDMPAnalyticSSFCTFilter(
  const std::shared_ptr<Limiter> limiter_, const DoFHandler<dim> & dof_handler_)
  : DMPSteadyStateFCTFilter<dim>(limiter_, dof_handler_),
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
void TransportDMPAnalyticSSFCTFilter<dim>::filter_antidiffusive_fluxes(
  const Vector<double> & solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs,
  const Vector<double> & cumulative_antidiffusion,
  SparseMatrix<double> & limiter_matrix,
  SparseMatrix<double> & antidiffusion_matrix)
{
  // compute DMP solution bounds
  compute_solution_bounds(solution, low_order_ss_matrix, ss_rhs);

  // compute analytic solution bounds
  analytic_bounds.update(solution, 0.0, 0.0);

  // widen solution bounds to include both DMP and analytic bounds
  this->solution_bounds.widen(analytic_bounds);

  // compute antidiffusion bounds Q- and Q+
  compute_antidiffusion_bounds(
    solution, low_order_ss_matrix, ss_rhs, cumulative_antidiffusion);

  // limit antidiffusion fluxes
  this->limiter->compute_limiter_matrix(
    antidiffusion_matrix, this->antidiffusion_bounds, limiter_matrix);
  this->limiter->apply_limiter_matrix(limiter_matrix, antidiffusion_matrix);
}