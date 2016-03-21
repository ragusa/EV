/**
 * \file DMPSteadyStateFCTFilter.cc
 * \brief Provides the function definitions for the DMPSteadyStateFCTFilter
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 */
template <int dim>
DMPSteadyStateFCTFilter<dim>::DMPSteadyStateFCTFilter(
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_)
  : SteadyStateFCTFilter<dim>(limiter_, dof_handler_, fe_)
{
}

/**
 * \brief Computes bounds to be imposed on the FCT solution, \f$W_i^-\f$ and
 *        \f$W_i^+\f$.
 *
 * The imposed bounds are computed as
 * \f[
 *   W_i^-(\mathrm{U}) = -\frac{1}{A^L_{i,j}}\sum\limits_{j\ne i}
 *     A^L_{i,j}U_{min,i} + \frac{b_i}{A^L_{i,i}} \,,
 * \f]
 * \f[
 *   W_i^+(\mathrm{U}) = -\frac{1}{A^L_{i,j}}\sum\limits_{j\ne i}
 *     A^L_{i,j}U_{max,i} + \frac{b_i}{A^L_{i,i}} \,,
 * \f]
 * where \f$U_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} U_j\f$ and
 * \f$U_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} U_j\f$.
 *
 * \param[in] solution  solution vector \f$\mathbf{U}\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 */
template <int dim>
void DMPSteadyStateFCTFilter<dim>::compute_solution_bounds(
  const Vector<double> & solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs)
{
  // create references
  Vector<double> & solution_min = this->solution_bounds.lower;
  Vector<double> & solution_max = this->solution_bounds.upper;

  // compute minimum and maximum values of solution
  this->compute_min_and_max_of_dof_vector(solution, solution_min, solution_max);

  // compute the upper and lower bounds for the FCT solution
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    double diagonal_term = 0.0;
    double off_diagonal_sum = 0.0;

    SparseMatrix<double>::const_iterator it = low_order_ss_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_end = low_order_ss_matrix.end(i);
    for (; it != it_end; ++it)
    {
      // get column index
      const unsigned int j = it->column();

      // get value
      const double ALij = it->value();

      // add nonzero entries to get the row sum
      if (j == i)
        diagonal_term = ALij;
      else
        off_diagonal_sum += ALij;
    }

    // compute bounds
    solution_max(i) = -off_diagonal_sum / diagonal_term * solution_max(i) +
      ss_rhs(i) / diagonal_term;
    solution_min(i) = -off_diagonal_sum / diagonal_term * solution_min(i) +
      ss_rhs(i) / diagonal_term;
  }
}

/**
 * \brief Computes the antidiffusion bounds \f$\mathbf{Q}^\pm\f$.
 *
 * \param[in] solution  solution \f$\mathbf{U}\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}\f$
 * \param[in] cumulative_antidiffusion  cumulative antidiffusion vector
 *            \f$\bar{\mathbf{p}}\f$
 */
template <int dim>
void DMPSteadyStateFCTFilter<dim>::compute_antidiffusion_bounds(
  const Vector<double> & solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs,
  const Vector<double> & cumulative_antidiffusion)
{
  // create references
  Vector<double> & Q_minus = this->antidiffusion_bounds.lower;
  Vector<double> & Q_plus = this->antidiffusion_bounds.upper;
  const Vector<double> & solution_min = this->solution_bounds.lower;
  const Vector<double> & solution_max = this->solution_bounds.upper;

  // initialize Q- and Q+
  Q_minus = 0;
  Q_minus.add(-1.0, ss_rhs);
  Q_minus.add(-1.0, cumulative_antidiffusion);
  Q_plus = 0;
  Q_plus.add(-1.0, ss_rhs);
  Q_plus.add(-1.0, cumulative_antidiffusion);

  // compute Q- and Q+
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    SparseMatrix<double>::const_iterator it = low_order_ss_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_end = low_order_ss_matrix.end(i);
    for (; it != it_end; ++it)
    {
      // get column index
      const unsigned int j = it->column();

      // get value
      const double ALij = it->value();

      // add to sums
      if (j == i) // diagonal element of A^L
      {
        Q_minus[i] += ALij * solution_min[i];
        Q_plus[i] += ALij * solution_max[i];
      }
      else // off-diagonal element of A^L
      {
        Q_minus[i] += ALij * solution[j];
        Q_plus[i] += ALij * solution[j];
      }
    }
  }
}
