/**
 * \file LinearSolver.cc
 * \brief Provides the function definitions for the LinearSolver class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
LinearSolver<dim>::LinearSolver(
  const LinearSolverType & linear_solver_type_,
  const ConstraintMatrix & constraints_,
  const DoFHandler<dim> & dof_handler_,
  std::shared_ptr<Function<dim>> dirichlet_function_,
  const unsigned int & n_components_)
  : linear_solver_type(linear_solver_type_),
    constraints(&constraints_),
    dof_handler(&dof_handler_),
    dirichlet_function(dirichlet_function_),
    have_dirichlet_bc(dirichlet_function != nullptr),
    n_components(n_components_)
{
}

/**
 * \brief Solves a linear system A*x = b, with no Dirichlet BC or constraints
 *        applied.
 *
 * \param[in]     A  system matrix
 * \param[in,out] x  system solution
 * \param[in]     b  system rhs
 */
template <int dim>
void LinearSolver<dim>::solve(const SparseMatrix<double> & A,
                              Vector<double> & x,
                              const Vector<double> & b) const
{
  // solve linear system
  switch (linear_solver_type)
  {
    case LinearSolverType::direct:
    {
      SparseDirectUMFPACK A_direct;
      A_direct.initialize(A);
      A_direct.vmult(x, b);
      break;
    }
    default:
    {
      throw ExcNotImplemented();
      break;
    }
  }
}

/**
 * \brief Solves a linear system A*x = b, with or without Dirichlet BC applied.
 *
 * \param[in] A  system matrix
 * \param[in,out] x  system solution
 * \param[in] b  system rhs
 * \param[in] dirichlet_bc_apply flag to apply Dirichlet BC
 * \param[in] t  time at which to evaluate Dirichlet BC function
 */
template <int dim>
void LinearSolver<dim>::solve_with_dirichlet(SparseMatrix<double> & A,
                                             Vector<double> & x,
                                             Vector<double> & b,
                                             const bool & dirichlet_bc_apply,
                                             const double & t)
{
  // apply Dirichlet BC if they apply
  if (dirichlet_bc_apply && have_dirichlet_bc)
    apply_dirichlet_bc(A, x, b, t);

  // set constrained DoFs to zero, see deal.II module on constrained degrees of
  // freedom - treatment of inhomogeneous constraints, approach 1; iterative
  // solvers might stop prematurely if this is not done
  constraints->set_zero(x);

  // solve linear system
  solve(A, x, b);

  // apply constraints
  constraints->distribute(x);
}

/**
 * \brief Applies Dirichlet boundary conditions to a linear system A*x = b.
 *
 * \param [in,out] A system matrix
 * \param [in,out] x system solution
 * \param [in,out] b system rhs
 * \param [in] t  time at which Dirichlet value function is to be evaluated
 */
template <int dim>
void LinearSolver<dim>::apply_dirichlet_bc(SparseMatrix<double> & A,
                                           Vector<double> & x,
                                           Vector<double> & b,
                                           const double & t)
{
  // create map of dofs to boundary values to be imposed
  std::map<unsigned int, double> boundary_values;

  // loop over components
  for (unsigned int component = 0; component < n_components; ++component)
  {
    // set time for dirichlet function
    dirichlet_function->set_time(t);

    // mask other components
    std::vector<bool> component_mask(n_components, false);
    component_mask[component] = true;

    // fill boundary_values with boundary values
    VectorTools::interpolate_boundary_values(
      *dof_handler, 0, *(dirichlet_function), boundary_values, component_mask);
  }

  // apply boundary values to system
  MatrixTools::apply_boundary_values(boundary_values, A, x, b);
}
