/**
 * \file DMPThetaFCTFilter.cc
 * \brief Provides the function definitions for the DMPThetaFCTFilter
 * class.
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
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 */
template <int dim>
DMPThetaFCTFilter<dim>::DMPThetaFCTFilter(
  const RunParameters & run_parameters_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const SparseMatrix<double> & lumped_mass_matrix_,
  const double & theta_,
  const std::map<unsigned int, double> & dirichlet_values_)
  : ThetaFCTFilter<dim>(run_parameters_,
                        limiter_,
                        dof_handler_,
                        fe_,
                        lumped_mass_matrix_,
                        theta_,
                        dirichlet_values_)
{
  // resize temporary vectors
  solution_min_new.reinit(this->n_dofs);
  solution_max_new.reinit(this->n_dofs);
  tmp_vector.reinit(this->n_dofs);
}

/**
 * \brief Computes bounds to be imposed on the FCT solution, \f$W_i^-\f$ and
 *        \f$W_i^+\f$.
 *
 * \warning This function assumes the conservation law flux is linear in that
 *          the low-order steady-state matrix \f$\mathbf{A}^L\f$ is assumed to
 *          not change with time.
 *
 * \param[in] new_solution  new solution vector \f$\mathbf{U}^{n+1}\f$
 * \param[in] old_solution  old solution vector \f$\mathbf{U}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs_new  new steady-state right hand side vector
 *            \f$\mathbf{b}^{n+1}\f$
 * \param[in] ss_rhs_old  old steady-state right hand side vector
 *            \f$\mathbf{b}^n\f$
 */
template <int dim>
void DMPThetaFCTFilter<dim>::compute_solution_bounds(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs_new,
  const Vector<double> & ss_rhs_old)
{
  // create references
  Vector<double> & solution_min = this->solution_bounds.lower;
  Vector<double> & solution_max = this->solution_bounds.upper;
  Vector<double> & solution_min_old = this->solution_bounds.lower;
  Vector<double> & solution_max_old = this->solution_bounds.upper;
  const SparseMatrix<double> & lumped_mass_matrix = *this->lumped_mass_matrix;

  // compute minimum and maximum values of solution
  this->compute_min_and_max_of_dof_vector(
    old_solution, solution_min_old, solution_max_old);
  this->compute_min_and_max_of_dof_vector(
    new_solution, solution_min_new, solution_max_new);

  // compute the upper and lower bounds for the FCT solution
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    // compute off diagonal row sum and diagonal term
    double off_diagonal_row_sum = 0.0;
    double diagonal_term = 0.0;

    SparseMatrix<double>::const_iterator it = low_order_ss_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_end = low_order_ss_matrix.end(i);
    for (; it != it_end; ++it)
    {
      // get column index
      const unsigned int j = it->column();

      // get value
      const double value = it->value();

      // add to sums
      if (j == i)
        diagonal_term = value;
      else
        off_diagonal_row_sum += value;
    }

    // compute full row sum
    const double row_sum = off_diagonal_row_sum + diagonal_term;

    // lumped mass matrix entry
    const double m_i = lumped_mass_matrix(i, i);

    // compute the max and min values for the maximum principle
    solution_max(i) =
      ((1.0 - (1.0 - this->theta) * dt / m_i * row_sum) * solution_max_old(i) -
       this->theta * dt / m_i * off_diagonal_row_sum * solution_max_new(i) +
       dt / m_i *
         ((1 - this->theta) * ss_rhs_old(i) + this->theta * ss_rhs_new(i))) /
      (1.0 + this->theta * dt / m_i * diagonal_term);
    solution_min(i) =
      ((1.0 - (1.0 - this->theta) * dt / m_i * row_sum) * solution_min_old(i) -
       this->theta * dt / m_i * off_diagonal_row_sum * solution_min_new(i) +
       dt / m_i *
         ((1 - this->theta) * ss_rhs_old(i) + this->theta * ss_rhs_new(i))) /
      (1.0 + this->theta * dt / m_i * diagonal_term);
  }
}

/**
 * \brief Computes the antidiffusion bounds \f$\mathbf{Q}^\pm\f$.
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
 */
template <int dim>
void DMPThetaFCTFilter<dim>::compute_antidiffusion_bounds(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs_new,
  const Vector<double> & ss_rhs_old,
  const Vector<double> & cumulative_antidiffusion)
{
  // create references
  Vector<double> & Q_minus = this->antidiffusion_bounds.lower;
  Vector<double> & Q_plus = this->antidiffusion_bounds.upper;
  const Vector<double> & solution_min = this->solution_bounds.lower;
  const Vector<double> & solution_max = this->solution_bounds.upper;
  Vector<double> & tmp = this->tmp_vector;
  const SparseMatrix<double> & lumped_mass_matrix = *this->lumped_mass_matrix;

  // compute Q+
  Q_plus = 0;
  lumped_mass_matrix.vmult(tmp, solution_max);
  Q_plus.add(1.0 / dt, tmp);
  lumped_mass_matrix.vmult(tmp, old_solution);
  Q_plus.add(-1.0 / dt, tmp);
  low_order_ss_matrix.vmult(tmp, old_solution);
  Q_plus.add(1.0 - this->theta, tmp);
  low_order_ss_matrix.vmult(tmp, new_solution);
  Q_plus.add(this->theta, tmp);
  Q_plus.add(-(1.0 - this->theta), ss_rhs_old);
  Q_plus.add(-this->theta, ss_rhs_new);
  Q_plus.add(-1.0, cumulative_antidiffusion);

  // compute Q-
  Q_minus = 0;
  lumped_mass_matrix.vmult(tmp, solution_min);
  Q_minus.add(1.0 / dt, tmp);
  lumped_mass_matrix.vmult(tmp, old_solution);
  Q_minus.add(-1.0 / dt, tmp);
  low_order_ss_matrix.vmult(tmp, old_solution);
  Q_minus.add(1.0 - this->theta, tmp);
  low_order_ss_matrix.vmult(tmp, new_solution);
  Q_minus.add(this->theta, tmp);
  Q_minus.add(-(1.0 - this->theta), ss_rhs_old);
  Q_minus.add(-this->theta, ss_rhs_new);
  Q_minus.add(-1.0, cumulative_antidiffusion);
}
