/**
 * \file ExplicitEulerFCTFilter.cc
 * \brief Provides the function definitions for the ExplicitEulerFCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 */
template <int dim>
ExplicitEulerFCTFilter<dim>::ExplicitEulerFCTFilter(
  const RunParameters & run_parameters_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const SparseMatrix<double> & lumped_mass_matrix_,
  const std::map<unsigned int, double> & dirichlet_values_)
  : FCTFilter<dim>(
      run_parameters_, limiter_, dof_handler_, fe_, dirichlet_values_),
    lumped_mass_matrix(&lumped_mass_matrix_)
{
  // resize temporary vector
  tmp_vector.reinit(this->n_dofs);
}

/**
 * \brief Computes the antidiffusion bounds \f$\mathbf{Q}^\pm\f$.
 *
 * \param[in] solution_bounds  solution bounds \f$\mathbf{W}^\pm\f$
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] inviscid_ss_flux  inviscid steady-state flux vector
 *   (entries are \f$(\mathbf{A}\mathbf{U}^n)_i\f$ for scalar case,
 *   \f$\sum\limits_j\mathbf{c}_{i,j}\cdot\mathrm{F}^n_j\f$ for systems case)
 * \param[in] low_order_diffusion_matrix  low-order diffusion matrix
 *            \f$\mathbf{D}^L\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}^n\f$
 */
template <int dim>
void ExplicitEulerFCTFilter<dim>::compute_antidiffusion_bounds(
  const DoFBounds<dim> & solution_bounds,
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> & inviscid_ss_flux,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const Vector<double> & ss_rhs)
{
  // create references
  Vector<double> & Q_minus = this->antidiffusion_bounds.lower;
  Vector<double> & Q_plus = this->antidiffusion_bounds.upper;
  const Vector<double> & solution_min = solution_bounds.lower;
  const Vector<double> & solution_max = solution_bounds.upper;
  Vector<double> & tmp = this->tmp_vector;
  const SparseMatrix<double> & lumped_mass_matrix = *this->lumped_mass_matrix;

  // start computing Q+
  Q_plus = 0;
  lumped_mass_matrix.vmult(tmp, old_solution);
  Q_plus.add(-1.0 / dt, tmp);
  Q_plus.add(1.0, inviscid_ss_flux);
  low_order_diffusion_matrix.vmult(tmp, old_solution);
  Q_plus.add(1.0, tmp);
  Q_plus.add(-1.0, ss_rhs);

  // copy current contents of Q+ as these components are identical
  Q_minus = Q_plus;

  // finish computing Q+ and Q-
  lumped_mass_matrix.vmult(tmp, solution_max);
  Q_plus.add(1.0 / dt, tmp);
  lumped_mass_matrix.vmult(tmp, solution_min);
  Q_minus.add(1.0 / dt, tmp);
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
void ExplicitEulerFCTFilter<dim>::filter_antidiffusive_fluxes(
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
  compute_antidiffusion_bounds(this->solution_bounds,
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
