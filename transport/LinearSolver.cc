/** \brief Constructor.
 */
template<int dim>
LinearSolver<dim>::LinearSolver(const unsigned int     &linear_solver_option,
                                const ConstraintMatrix &constraints,
                                const DoFHandler<dim>  &dof_handler,
                                const Function<dim>    &dirichlet_value_function) :
   linear_solver_option(linear_solver_option),
   constraints(constraints),
   dof_handler(&dof_handler),
   dirichlet_value_function(&dirichlet_value_function)
{
}

/** \brief Destructor.
 */
template<int dim>
LinearSolver<dim>::~LinearSolver()
{
}

/** \brief Solve a linear system A*x = b.
 *  \param [in] A system matrix
 *  \param [in] b system rhs
 *  \param [in,out] x system solution
 */
template<int dim>
void LinearSolver<dim>::solve(SparseMatrix<double> &A,
                              Vector<double>       &b,
                              Vector<double>       &x) const
{
   // apply Dirichlet BC
   apply_Dirichlet_BC(A,b,x);

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
   constraints.distribute(x);
}

/** \brief Applies Dirichlet boundary conditions to a linear system A*x = b.
 *  \param [in,out] A system matrix
 *  \param [in,out] b system rhs
 *  \param [in,out] x system solution
 */
template<int dim>
void LinearSolver<dim>::apply_Dirichlet_BC(SparseMatrix<double> &A,
                                           Vector<double>       &b,
                                           Vector<double>       &x) const
{
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
