/**
 * \file SteadyStateFCTFilter.cc
 * \brief Provides the function definitions for the SteadyStateFCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 */
template <int dim>
SteadyStateFCTFilter<dim>::SteadyStateFCTFilter(
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_)
  : FCTFilter<dim>(limiter_, dof_handler_)
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
void SteadyStateFCTFilter<dim>::filter_antidiffusive_fluxes(
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
  compute_antidiffusion_bounds(
    solution, low_order_ss_matrix, ss_rhs, cumulative_antidiffusion);

  // limit antidiffusion fluxes
  this->limiter->compute_limiter_matrix(
    antidiffusion_matrix, this->antidiffusion_bounds, limiter_matrix);
  this->limiter->apply_limiter_matrix(limiter_matrix, antidiffusion_matrix);
}
