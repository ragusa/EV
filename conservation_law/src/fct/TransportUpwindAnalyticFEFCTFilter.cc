/**
 * \file TransportUpwindAnalyticFEFCTFilter.cc
 * \brief Provides the function definitions for the
 *        TransportUpwindAnalyticFEFCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] problem_parameters_  problem parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] limiter_  limiter
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 */
template <int dim>
TransportUpwindAnalyticFEFCTFilter<dim>::TransportUpwindAnalyticFEFCTFilter(
  const RunParameters & run_parameters_,
  TransportProblemParameters<dim> & problem_parameters_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const SparseMatrix<double> & lumped_mass_matrix_,
  const std::map<unsigned int, double> & dirichlet_values_)
  : ExplicitEulerFCTFilter<dim>(run_parameters_,
                                limiter_,
                                dof_handler_,
                                fe_,
                                lumped_mass_matrix_,
                                dirichlet_values_),
    analytic_bounds(problem_parameters_,
                    dof_handler_,
                    fe_,
                    run_parameters_.upwind_bounds_sampling_points),
    speed(problem_parameters_.transport_speed)
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
 * \param[in] t_old  old time
 * \param[inout] limiter_matrix  limiter matrix \f$\mathbf{L}\f$
 * \param[inout] antidiffusion_matrix  antidiffusion matrix \f$\mathbf{P}\f$
 */
template <int dim>
void TransportUpwindAnalyticFEFCTFilter<dim>::filter_antidiffusive_fluxes(
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> & inviscid_ss_flux,
  const Vector<double> & ss_reaction,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const Vector<double> & ss_rhs,
  const double & t_old,
  SparseMatrix<double> & limiter_matrix,
  SparseMatrix<double> & antidiffusion_matrix)
{
  // compute solution bounds W- and W+
  compute_solution_bounds(old_solution, dt, ss_reaction, ss_rhs, t_old);

  // compute antidiffusion bounds Q- and Q+
  this->compute_antidiffusion_bounds(analytic_bounds,
                                     old_solution,
                                     dt,
                                     inviscid_ss_flux,
                                     low_order_diffusion_matrix,
                                     ss_rhs);

  // enforce antidiffusion bounds signs if requested
  if (this->do_enforce_antidiffusion_bounds_signs)
    this->enforce_antidiffusion_bounds_signs();

  // check signs of antidiffusion bounds
  this->check_antidiffusion_bounds_signs();

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
bool TransportUpwindAnalyticFEFCTFilter<dim>::check_bounds(
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
Vector<double> TransportUpwindAnalyticFEFCTFilter<dim>::get_lower_solution_bound()
  const
{
  return analytic_bounds.lower;
}

/**
 * \brief Gets the upper solution bound.
 *
 * \return upper solution bound vector
 */
template <int dim>
Vector<double> TransportUpwindAnalyticFEFCTFilter<dim>::get_upper_solution_bound()
  const
{
  return analytic_bounds.upper;
}

/**
 * \brief Computes bounds to be imposed on the FCT solution, \f$W_i^-\f$ and
 *        \f$W_i^+\f$.
 *
 * \param[in] old_solution  old solution vector \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] ss_reaction  steady-state reaction vector \f$\mathbf{\sigma}\f$,
 *            where \f$\sigma_i = \int\limits_{S_i}
 *            \varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *            = \sum\limits_j\int\limits_{S_{i,j}}
 *            \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 * \param[in] t_old  old time
 */
template <int dim>
void TransportUpwindAnalyticFEFCTFilter<dim>::compute_solution_bounds(
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> &,
  const Vector<double> &,
  const double & t_old)
{
  // compute analytic solution bounds
  analytic_bounds.update(old_solution, speed * dt, t_old);
}
