/**
 * \file DMPExplicitEulerFCTFilter.cc
 * \brief Provides the function definitions for the DMPExplicitEulerFCTFilter
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 */
template <int dim>
DMPExplicitEulerFCTFilter<dim>::DMPExplicitEulerFCTFilter(
  const std::shared_ptr<Limiter> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const SparseMatrix<double> & lumped_mass_matrix_)
  : ExplicitEulerFCTFilter<dim>(limiter_, dof_handler_, lumped_mass_matrix_)
{
  // resize temporary vector
  tmp_vector.reinit(this->n_dofs);
}

/**
 * \brief Computes bounds to be imposed on the FCT solution, \f$W_i^-\f$ and
 *        \f$W_i^+\f$.
 *
 * The imposed bounds are computed as
 * \f[
 *   W_i^-(\mathrm{U}^n) = U_{min,i}^n\left(1
 *     - \frac{\Delta t}{M^L_{i,i}}\sigma_i\right)
 *     + \frac{\Delta t}{M^L_{i,i}}b_i^n \,,
 * \f]
 * \f[
 *   W_i^+(\mathrm{U}^n) = U_{max,i}^n\left(1
 *     - \frac{\Delta t}{M^L_{i,i}}\sigma_i\right)
 *     + \frac{\Delta t}{M^L_{i,i}}b_i^n \,,
 * \f]
 * where \f$U_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} U_j\f$ and
 * \f$U_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} U_j\f$.
 *
 * \param[in] old_solution  old solution vector \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] ss_reaction  steady-state reaction vector \f$\mathbf{\sigma}\f$,
 *            where \f$\sigma_i = \int\limits_{S_i}
 *            \varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *            = \sum\limits_j\int\limits_{S_{i,j}}
 *            \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 */
template <int dim>
void DMPExplicitEulerFCTFilter<dim>::compute_solution_bounds(
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> & ss_reaction,
  const Vector<double> & ss_rhs)
{
  // create references
  Vector<double> & solution_min = this->solution_bounds.lower;
  Vector<double> & solution_max = this->solution_bounds.upper;
  const SparseMatrix<double> & lumped_mass_matrix = *this->lumped_mass_matrix;

  // compute minimum and maximum values of solution
  this->compute_min_and_max_of_dof_vector(
    old_solution, solution_min, solution_max);

  /*
    // consider star states in bounds if specified
    if (use_star_states_in_fct_bounds)
    {
      for (unsigned int i = 0; i < n_dofs; ++i)
      {
        // iterate over sparsity pattern to get degree of freedom indices
        // in the support of i
        SparsityPattern::iterator it = sparsity->begin(i);
        SparsityPattern::iterator it_end = sparsity->end(i);
        for (; it != it_end; ++it)
        {
          // get column index
          const unsigned int j = it->column();

          if (j != i)
          {
            // determine if j corresponds to the same component as i
            if (i % n_components == j % n_components)
            {
              // get star state value associated with i and j
              const double star_state_value = star_state->get_star_state(i, j);
              solution_max(i) = std::max(solution_max(i), star_state_value);
              solution_min(i) = std::min(solution_min(i), star_state_value);
            }
          }
        }
      }
    }
  */

  // At this point, the min/max values of the old solution in the support
  // of test function i are stored in solution_min(i) and solution_max(i).
  // Now these values are multiplied by (1-dt/m(i))*sigma(i) and
  // added to dt/m(i)*b(i).

  // compute the upper and lower bounds for the FCT solution
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    solution_max(i) =
      solution_max(i) * (1.0 - dt / lumped_mass_matrix(i, i) * ss_reaction(i)) +
      dt / lumped_mass_matrix(i, i) * ss_rhs(i);

    solution_min(i) =
      solution_min(i) * (1.0 - dt / lumped_mass_matrix(i, i) * ss_reaction(i)) +
      dt / lumped_mass_matrix(i, i) * ss_rhs(i);
  }
}

/**
 * \brief Computes the antidiffusion bounds \f$\mathbf{Q}^\pm\f$.
 *
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
void DMPExplicitEulerFCTFilter<dim>::compute_antidiffusion_bounds(
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> & inviscid_ss_flux,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const Vector<double> & ss_rhs)
{
  // create references
  Vector<double> & Q_minus = this->antidiffusion_bounds.lower;
  Vector<double> & Q_plus = this->antidiffusion_bounds.upper;
  const Vector<double> & solution_min = this->solution_bounds.lower;
  const Vector<double> & solution_max = this->solution_bounds.upper;
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
