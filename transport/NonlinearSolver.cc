/**
 * \brief Constructor.
 *
 * \param[in] parameters_  input parameters
 * \param[in] constraints_  constraints for linear system
 * \param[in] dof_handler_  dof handler
 * \param[in] dirichlet_value_function_  function for the Dirichlet BC to be
 *   imposed on the linear system
 */
template<int dim>
NonlinearSolver<dim>::NonlinearSolver(
  const TransportParameters<dim> & parameters_,
  const ConstraintMatrix         & constraints_,
  const DoFHandler<dim>          & dof_handler_,
  Function<dim>                  & dirichlet_value_function_) :
    relaxation_factor(parameters_.relaxation_factor),
    nonlinear_tolerance(parameters_.nonlinear_tolerance),
    iteration_max(parameters_.nonlinear_max_iteration),
    constraints(constraints_),
    dof_handler(dof_handler_),
    dirichlet_value_function(dirichlet_value_function_),
    linear_solver(parameters_.linear_solver_option,
      constraints,
      dof_handler,
      dirichlet_value_function),
    iteration_number(0)
{
  // get number of degrees of freedom
  n_dofs = dof_handler_.n_dofs();

  // reinitialize vectors
  solution_change.reinit(n_dofs);
  residual.reinit(n_dofs);
}

/**
 * \brief Initializes the solution to guess and sets iteration counter to zero.
 *
 * \param[in] solution_guess  initial guess for the solution
 */
template<int dim>
void NonlinearSolver<dim>::initialize(Vector<double> & solution_guess)
{
  // get pointer to solution vector, which current has a guess for the solution
  solution = &solution_guess;

  // set the interation number to zero
  iteration_number = 0;

  std::cout << "  Nonlinear solve:" << std::endl;
}

/**
 * \brief Checks convergence of nonlinear system.
 *
 * The L-2 norm of the linear residual is used to determine convergence.
 *
 * \param[in] residual  residual for linear system
 * \return boolean convergence flag
 */
template<int dim>
bool NonlinearSolver<dim>::check_convergence(const Vector<double> & residual)
{
  // increment iteration number
  ++iteration_number;

  // compute L-2 norm of linear residual
  double nonlinear_err = residual.l2_norm();

  // print error
  std::cout << "  (" << iteration_number << ") err = " <<
    nonlinear_err << std::endl;

  // determine if error is within the nonlinear tolerance
  if (nonlinear_err < nonlinear_tolerance)
  {
    return true;
  }
  else
  {
    // check if max iteration has been reached
    AssertThrow(iteration_number < iteration_max,
      ExcMaxIterationReached(iteration_max));

    return false;
  }
}

/**
 * \brief Computes the linear residual, checks convergence, and if not
 *        converged, updates the solution iterate.
 *
 * \param[in] A system matrix
 * \param[in] b system right-hand-side
 *
 * \return convergence flag
 */
template<int dim>
bool NonlinearSolver<dim>::update(
  const SparseMatrix<double> & A,
  const Vector<double> & b)
{
  // compute linear residual: r^(l) = b^(l) - A^(l)*U^(l)
  A.vmult(residual, *solution);
  residual *= -1.0;
  residual.add(1.0, b);

  // check convergence
  bool converged = check_convergence(residual);

  // solve for solution update dU: A^(l)*dU = r^(l)
  linear_solver.solve(A, residual, solution_change);
  //linear_solver.solve(A, solution_change, residual);

  // update solution: U^(l+1) = U^(l) + relax*dU
  solution->add(relaxation_factor, solution_change);

  // distribute constraints
  constraints.distribute(*solution);

  // return convergence flag
  return converged;
}

