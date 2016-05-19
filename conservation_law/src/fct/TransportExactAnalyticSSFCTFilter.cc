/**
 * \file TransportExactAnalyticSSFCTFilter.cc
 * \brief Provides the function definitions for the
 *        TransportExactAnalyticSSFCTFilter class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] problem_parameters_  problem parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] cell_quadrature_  cell quadrature
 * \param[in] limiter_  limiter
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 */
template <int dim>
TransportExactAnalyticSSFCTFilter<dim>::TransportExactAnalyticSSFCTFilter(
  const RunParameters & run_parameters_,
  TransportProblemParameters<dim> & problem_parameters_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const QGauss<dim> & cell_quadrature_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const std::map<unsigned int, double> & dirichlet_values_)
  : TransportAnalyticSSFCTFilter<dim>(run_parameters_,
                                      problem_parameters_,
                                      dof_handler_,
                                      fe_,
                                      cell_quadrature_,
                                      limiter_,
                                      dirichlet_values_)
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
 * \brief Computes bounds to be imposed on the FCT solution, \f$W_i^-\f$ and
 *        \f$W_i^+\f$.
 *
 * \param[in] solution  solution vector \f$\mathbf{U}\f$
 * \param[in] low_order_ss_matrix  low-order steady-state matrix
 *            \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 */
template <int dim>
void TransportExactAnalyticSSFCTFilter<dim>::compute_solution_bounds(
  const Vector<double> &, const SparseMatrix<double> &, const Vector<double> &)
{
  // compute analytic solution bounds
  this->analytic_bounds.update(exact_solution_values, 0.0, 0.0);
}
