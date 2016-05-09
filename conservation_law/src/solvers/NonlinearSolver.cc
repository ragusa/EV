/**
 * \brief Constructor.
 *
 * \param[in] parameters_     input parameters
 * \param[in] linear_solver_  linear solver
 * \param[in] constraints_    constraints for linear system
 */
template <int dim>
NonlinearSolver<dim>::NonlinearSolver(const RunParameters & parameters_,
                                      const LinearSolver<dim> & linear_solver_,
                                      const ConstraintMatrix & constraints_)
  : relaxation_factor(parameters_.relaxation_factor),
    nonlinear_tolerance(parameters_.nonlinear_tolerance),
    iteration_max(parameters_.nonlinear_max_iterations),
    linear_solver(&linear_solver_),
    constraints(constraints_),
    iteration_number(0),
    n_dofs(0),
    cout1(std::cout, parameters_.verbosity_level >= 1)
{
}

/**
 * \brief Reinitializes vectors with new system size.
 *
 * \param[in] n_dofs  number of degrees of freedom
 */
template <int dim>
void NonlinearSolver<dim>::reinit(const unsigned int & n_dofs)
{
  // reinitialize vectors
  solution_change.reinit(n_dofs);
  residual.reinit(n_dofs);
}

/**
 * \brief Initializes the solution to guess and sets iteration counter to zero.
 *
 * \param[in] solution_guess  initial guess for the solution
 */
template <int dim>
void NonlinearSolver<dim>::initialize(Vector<double> & solution_guess)
{
  // get pointer to solution vector, which current has a guess for the solution
  solution = &solution_guess;

  // set the interation number to zero
  iteration_number = 0;

  // update system size if necessary
  if (solution->size() != n_dofs)
  {
    // update current system size
    n_dofs = solution->size();

    // reinitialize vectors
    reinit(n_dofs);
  }

  cout1 << "  Nonlinear solve:" << std::endl;
}

/**
 * \brief Checks convergence of nonlinear system.
 *
 * The L-2 norm of the linear residual is used to determine convergence.
 *
 * \param[in] residual  residual for linear system
 * \return boolean convergence flag
 */
template <int dim>
bool NonlinearSolver<dim>::check_convergence(const Vector<double> & residual)
{
  // increment iteration number
  ++iteration_number;

  // compute L-2 norm of linear residual
  double nonlinear_err = residual.l2_norm();

  // print error
  cout1 << "  (" << iteration_number << ") err = " << nonlinear_err << std::endl;

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
template <int dim>
bool NonlinearSolver<dim>::update(const SparseMatrix<double> & A,
                                  const Vector<double> & b)
{
  // compute linear residual: r^(l) = b^(l) - A^(l)*U^(l)
  A.vmult(residual, *solution);
  residual *= -1.0;
  residual.add(1.0, b);

  // check convergence
  bool converged = check_convergence(residual);
  if (converged)
    return converged;

  // solve for solution update dU: A^(l)*dU = r^(l)
  linear_solver->solve(A, solution_change, residual);

  // update solution: U^(l+1) = U^(l) + relax*dU
  solution->add(relaxation_factor, solution_change);

  // distribute constraints
  constraints.distribute(*solution);

  // return convergence flag
  return converged;
}

/**
 * \brief Returns the number of nonlinear iterations since the last
 *        \c initialize() call
 *
 * \return  number of nonlinear iterations since the last \c initialize() call
 */
template <int dim>
unsigned int NonlinearSolver<dim>::get_number_of_iterations() const
{
  return iteration_number;
}
