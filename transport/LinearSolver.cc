/** \brief Constructor.
 */
template<int dim>
LinearSolver<dim>::LinearSolver(const unsigned int     & linear_solver_option,
                                const ConstraintMatrix & constraints,
                                const DoFHandler<dim>  & dof_handler,
                                Function<dim>          & dirichlet_value_function) :
   linear_solver_option(linear_solver_option),
   constraints(& constraints),
   dof_handler(& dof_handler),
   dirichlet_value_function(& dirichlet_value_function)
{
}

/** \brief Destructor.
 */
template<int dim>
LinearSolver<dim>::~LinearSolver()
{
}

/**
 * Solves a linear system A*x = b, with or without Dirichlet BC applied.
 *
 * @param[in] A  system matrix
 * @param[in] b  system rhs
 * @param[in,out] x  system solution
 * @param[in] apply_dirichlet_bc  option to apply Dirichlet BC
 * @param[in] t  time at which to evaluate Dirichlet BC function
 */
template<int dim>
void LinearSolver<dim>::solve(
  SparseMatrix<double> & A,
  Vector<double>       & b,
  Vector<double>       & x,
  const bool           & apply_dirichlet_bc,
  const double         & t)
{
   // apply Dirichlet BC if they apply
   if (apply_dirichlet_bc)
     apply_Dirichlet_BC(A,b,x,t);

   // solve linear system
   switch (linear_solver_option) {
      case 1: {
         SparseDirectUMFPACK A_direct;
         A_direct.initialize(A);
         A_direct.vmult(x, b);
         break;
      } default: {
         Assert(false,ExcNotImplemented());
         break;
      }
   }

   // apply constraints
   constraints->distribute(x);
}

/**
 * Applies Dirichlet boundary conditions to a linear system A*x = b.
 *
 * @param [in,out] A system matrix
 * @param [in,out] b system rhs
 * @param [in,out] x system solution
 * @param [in] t  time at which Dirichlet value function is to be evaluated
 */
template<int dim>
void LinearSolver<dim>::apply_Dirichlet_BC(SparseMatrix<double> &A,
                                           Vector<double>       &b,
                                           Vector<double>       &x,
                                           const double          t)
{
   // set time for dirichlet function
   dirichlet_value_function->set_time(t);

   // create map of dofs to boundary values to be imposed
   std::map<unsigned int, double> boundary_values;
   VectorTools::interpolate_boundary_values(*dof_handler,
                                            1,
                                            *dirichlet_value_function,
                                            boundary_values);
   // apply boundary values to system
   MatrixTools::apply_boundary_values(boundary_values,
                                      A,
                                      x,
                                      b);
}
