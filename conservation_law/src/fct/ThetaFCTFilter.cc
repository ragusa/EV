/**
 * \file ThetaFCTFilter.cc
 * \brief Provides the function definitions for the ThetaFCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 * \param[in] theta_  theta parameter \f$\theta\f$
 */
template <int dim>
ThetaFCTFilter<dim>::ThetaFCTFilter(
  const RunParameters & run_parameters_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const SparseMatrix<double> & lumped_mass_matrix_,
  const double & theta_)
  : FCTFilter<dim>(run_parameters_, limiter_, dof_handler_, fe_),
    lumped_mass_matrix(&lumped_mass_matrix_),
    theta(theta_)
{
}

/**
 * \brief Filters antidiffusive fluxes.
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
void ThetaFCTFilter<dim>::filter_antidiffusive_fluxes(
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
  // compute solution bounds W- and W+
  compute_solution_bounds(
    new_solution, old_solution, dt, low_order_ss_matrix, ss_rhs_new, ss_rhs_old);

  // compute antidiffusion bounds Q- and Q+
  compute_antidiffusion_bounds(new_solution,
                               old_solution,
                               dt,
                               low_order_ss_matrix,
                               ss_rhs_new,
                               ss_rhs_old,
                               cumulative_antidiffusion);

  // enforce antidiffusion bounds signs if requested
  if (this->do_enforce_antidiffusion_bounds_signs)
    this->enforce_antidiffusion_bounds_signs();

  // limit antidiffusion fluxes
  this->limiter->compute_limiter_matrix(
    antidiffusion_matrix, this->antidiffusion_bounds, limiter_matrix);
  this->limiter->apply_limiter_matrix(limiter_matrix, antidiffusion_matrix);
}
