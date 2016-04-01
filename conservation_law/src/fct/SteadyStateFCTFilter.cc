/**
 * \file SteadyStateFCTFilter.cc
 * \brief Provides the function definitions for the SteadyStateFCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 */
template <int dim>
SteadyStateFCTFilter<dim>::SteadyStateFCTFilter(
  const RunParameters & run_parameters_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_)
  : FCTFilter<dim>(run_parameters_, limiter_, dof_handler_, fe_)
{
}

/**
 * \brief Computes the antidiffusion bounds \f$\mathbf{Q}^\pm\f$.
 *
 * \param[in] solution_bounds  solution bounds \f$\mathbf{W}^\pm\f$
 * \param[in] solution  solution \f$\mathbf{U}\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}\f$
 * \param[in] cumulative_antidiffusion  cumulative antidiffusion vector
 *            \f$\bar{\mathbf{p}}\f$
 */
template <int dim>
void SteadyStateFCTFilter<dim>::compute_antidiffusion_bounds(
  const DoFBounds<dim> & solution_bounds,
  const Vector<double> & solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs,
  const Vector<double> & cumulative_antidiffusion)
{
  // create references
  Vector<double> & Q_minus = this->antidiffusion_bounds.lower;
  Vector<double> & Q_plus = this->antidiffusion_bounds.upper;
  const Vector<double> & solution_min = solution_bounds.lower;
  const Vector<double> & solution_max = solution_bounds.upper;

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

  // enforce antidiffusion bounds signs if requested
  if (this->do_enforce_antidiffusion_bounds_signs)
    this->enforce_antidiffusion_bounds_signs();
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
  compute_antidiffusion_bounds(this->solution_bounds,
                               solution,
                               low_order_ss_matrix,
                               ss_rhs,
                               cumulative_antidiffusion);

  // limit antidiffusion fluxes
  this->limiter->compute_limiter_matrix(
    antidiffusion_matrix, this->antidiffusion_bounds, limiter_matrix);
  this->limiter->apply_limiter_matrix(limiter_matrix, antidiffusion_matrix);
}
