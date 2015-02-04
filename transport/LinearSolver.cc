/** \brief Constructor.
 */
LinearSolver::LinearSolver(const unsigned int     &linear_solver_option,
                           const ConstraintMatrix &constraints) :
   linear_solver_option(linear_solver_option),
   constraints(constraints)
{
}

/** \brief Destructor.
 */
LinearSolver::~LinearSolver()
{
}

/** \brief Solve a linear system A*x = b.
 */
void LinearSolver::solve(const SparseMatrix<double> &A,
                         const Vector<double>       &b,
                         Vector<double>             &x) const
{
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
