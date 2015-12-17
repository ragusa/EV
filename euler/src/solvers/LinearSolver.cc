/**
 * \brief Constructor.
 */
template <int dim>
LinearSolver<dim>::LinearSolver(
  const LinearSolverType & linear_solver_type_,
  const ConstraintMatrix & constraints_,
  const DoFHandler<dim> & dof_handler_,
  std::vector<FunctionParser<dim> *> dirichlet_functions_)
  : linear_solver_type(linear_solver_type_),
    constraints(&constraints_),
    dof_handler(&dof_handler_),
    dirichlet_functions(dirichlet_functions_),
    n_dirichlet_boundaries(dirichlet_functions_.size()),
    have_dirichlet_bc(n_dirichlet_boundaries > 0),
    n_components(dof_handler_.get_fe().n_components())
{
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
void LinearSolver<dim>::solve(SparseMatrix<double> & A,
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
      Assert(false, ExcNotImplemented());
      break;
    }
  }

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

  // loop over boundary IDs
  for (unsigned int boundary = 0; boundary < n_dirichlet_boundaries; ++boundary)
    // loop over components
    for (unsigned int component = 0; component < n_components; ++component)
    {
      // set time for dirichlet function
      dirichlet_functions[boundary]->set_time(t);

      // mask other components
      std::vector<bool> component_mask(n_components, false);
      component_mask[component] = true;

      // fill boundary_values with boundary values
      VectorTools::interpolate_boundary_values(*dof_handler,
                                               boundary,
                                               *(dirichlet_functions[boundary]),
                                               boundary_values,
                                               component_mask);
    }

  // apply boundary values to system
  MatrixTools::apply_boundary_values(boundary_values, A, x, b);
}
