/**
 * Constructor.
 *
 * @param[in] parameters_  input parameters
 * @param[in] constraints_  constraints for linear system
 * @param[in] dof_handler_  dof handler
 * @param[in] dirichlet_value_function_  function for the Dirichlet BC to be
 *   imposed on the linear system
 */
template<int dim>
NonlinearSolver<dim>::NonlinearSolver(
  const TransportParameters<dim> & parameters_,
  const ConstraintMatrix         & constraints_,
  const DoFHandler<dim>          & dof_handler_,
  Function<dim>                  & dirichlet_value_function_) :
    nonlinear_tolerance(parameters_.nonlinear_tolerance),
    iteration_max(parameters_.nonlinear_iteration_max),
    constraints(constraints_),
    dof_handler(dof_handler_),
    dirichlet_value_function(dirichlet_value_function_),
    linear_solver(parameters_.linear_solver_option,
      constraints,
      dof_handler,
      dirichlet_value_function),
    iteration_number(0)
{
}

/**
 * Resets the nonlinear solver.
 *
 * @param[in] solution_guess  initial guess for the solution
 */
template<int dim>
void NonlinearSolver<dim>::reset(const Vector<double> & solution_guess)
{
  // initialize the guess for the solution
  solution = solution_guess;

  // set the interation number to zero
  iteration_number = 0;
}

/**
 * Checks convergence of nonlinear system.
 *
 * @param[in] new_solution  new solution
 * @return boolean convergence flag
 */
template<int dim>
bool NonlinearSolver<dim>::checkConvergence(const Vector<double> & new_solution)
{
  // increment iteration number
  ++iteration_number;

  // compute L-2 norm of difference between new solution and previous solution
  Vector<double> & difference = solution;
  difference.add(-1.0, new_solution);
  double nonlinear_err = difference.l2_norm();

  // print error
  std::cout << "Iter " << iteration_number << ": err = " <<
    nonlinear_err << std::endl;

  // determine if error is within the nonlinear tolerance
  if (nonlinear_err < nonlinear_tolerance)
    return true;
  else
  {
    // check if max iteration has been reached
    if (iteration_number == iteration_max)
      ExcMaxIterationReached(iteration_max);

    // keep new solution as next iteration's previous iterate
    solution = new_solution;

    return false;
  }
}

/**
 * Gets the solution.
 *
 * @return solution vector
 */
template<int dim>
Vector<double> NonlinearSolver<dim>::getSolution() const
{
  return solution;
}
