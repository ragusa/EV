/**
 * \file TransportExactDMPSSFCTFilter.cc
 * \brief Provides the function definitions for the TransportExactDMPSSFCTFilter
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 * \param[in] problem_parameters_  problem parameters
 */
template <int dim>
TransportExactDMPSSFCTFilter<dim>::TransportExactDMPSSFCTFilter(
  const RunParameters & run_parameters_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const std::map<unsigned int, double> & dirichlet_values_,
  TransportProblemParameters<dim> & problem_parameters_)
  : SteadyStateFCTFilter<dim>(
      run_parameters_, limiter_, dof_handler_, fe_, dirichlet_values_),
    dmp_filter(run_parameters_, limiter_, dof_handler_, fe_, dirichlet_values_)
{
  // resize exact solution values vector
  exact_solution_values.reinit(dof_handler_.n_dofs());

  // assert that problem has an exact solution
  AssertThrow(problem_parameters_.has_exact_solution, ExcNoExactSolution());

  // compute exact solution values at each node
  VectorTools::interpolate(dof_handler_,
                           *(problem_parameters_.exact_solution_function),
                           exact_solution_values);
}

/**
 * \brief Evaluates the DMP solution bounds \f$W_i^-\f$ and \f$W_i^+\f$
 *        with the exact solution.
 *
 * The DMP solution bounds are computed as
 * \f[
 *   W_i^-(\mathrm{U}) = -\frac{1}{A^L_{i,j}}\sum\limits_{j\ne i}
 *     A^L_{i,j}U_{min,i} + \frac{b_i}{A^L_{i,i}} \,,
 * \f]
 * \f[
 *   W_i^+(\mathrm{U}) = -\frac{1}{A^L_{i,j}}\sum\limits_{j\ne i}
 *     A^L_{i,j}U_{max,i} + \frac{b_i}{A^L_{i,i}} \,,
 * \f]
 * where instead of
 * \f[
 *   U_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} U_j
 * \f]
 * \f[
 *   U_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} U_j
 * \f]
 * the exact solution is used in place of solution values:
 * \f[
 *   U_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} u^{exact}(\mathbf{x}_j)
 * \f]
 * \f[
 *   U_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} u^{exact}(\mathbf{x}_j)
 * \f]
 *
 * \param[in] solution  solution vector \f$\mathbf{U}\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 */
template <int dim>
void TransportExactDMPSSFCTFilter<dim>::compute_solution_bounds(
  const Vector<double> &,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs)
{
  // evaluate the DMP solution bounds with the exact solution values
  dmp_filter.compute_solution_bounds(
    exact_solution_values, low_order_ss_matrix, ss_rhs);

  // copy bounds from DMP filter object
  this->solution_bounds.lower = dmp_filter.get_lower_solution_bound();
  this->solution_bounds.upper = dmp_filter.get_upper_solution_bound();
}
