/**
 * Constructor.
 *
 * @param[in] parameters input parameters
 */
template<int dim>
NonlinearSolver<dim>::NonlinearSolver(
  const TransportParameters<dim> & parameters_) :
    nonlinear_tolerance(parameters_.nonlinear_tolerance)
{
}

/**
 * Ends an iteration.
 */
template<int dim>
void NonlinearSolver<dim>::endIteration()
{
  if (iteration_number
}

/**
 * Checks convergence of nonlinear system.
 *
 * @param[in] new_solution new iterate of solution vector
 * @return boolean convergence flag
 */
template<int dim>
bool NonlinearSolver<dim>::checkConvergence(const Vector<double> & new_solution)
{
  // compute L-2 norm of difference between new solution and previous solution
  Vector<double> & difference = previous_solution;
  difference.add(-1.0, new_solution);
  double nonlinear_err = difference.l2_norm();

  // keep new solution as next iteration's previous iterate
  previous_solution = new_solution;

  // determine if error is within the nonlinear tolerance
  if (nonlinear_err < nonlinear_tolerance)
    return true;
  else
    return false;
}
