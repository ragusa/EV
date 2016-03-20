/**
 * \file ExplicitEulerFCTFilter.cc
 * \brief Provides the function definitions for the ExplicitEulerFCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 */
template <int dim>
ExplicitEulerFCTFilter<dim>::ExplicitEulerFCTFilter(
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const SparseMatrix<double> & lumped_mass_matrix_)
  : FCTFilter<dim>(limiter_, dof_handler_),
    lumped_mass_matrix(&lumped_mass_matrix_)
{
}

/**
 * \brief Filters antidiffusive fluxes.
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
void ExplicitEulerFCTFilter<dim>::filter_antidiffusive_fluxes(
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> & inviscid_ss_flux,
  const Vector<double> & ss_reaction,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const Vector<double> & ss_rhs,
  SparseMatrix<double> & limiter_matrix,
  SparseMatrix<double> & antidiffusion_matrix)
{
  // compute solution bounds W- and W+
  compute_solution_bounds(old_solution, dt, ss_reaction, ss_rhs);

  // compute antidiffusion bounds Q- and Q+
  compute_antidiffusion_bounds(
    old_solution, dt, inviscid_ss_flux, low_order_diffusion_matrix, ss_rhs);

  // limit antidiffusion fluxes
  this->limiter->compute_limiter_matrix(
    antidiffusion_matrix, this->antidiffusion_bounds, limiter_matrix);
  this->limiter->apply_limiter_matrix(limiter_matrix, antidiffusion_matrix);
}
