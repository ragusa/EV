#ifndef LinearSolver_cc
#define LinearSolver_cc

#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

using namespace dealii;

/**
 * Class for solving linear systems.
 */
template<int dim>
class LinearSolver {
   public:
      LinearSolver(const unsigned int     & linear_solver_option,
                   const ConstraintMatrix & constraints,
                   const DoFHandler<dim>  & dof_handler,
                   Function<dim>          & dirichlet_value_function);
      ~LinearSolver();

      void solve(const SparseMatrix<double> & A,
                 const Vector<double>       & b,
                 Vector<double>             & x);

      void solve(SparseMatrix<double> & A,
                 Vector<double>       & b,
                 Vector<double>       & x,
                 const bool           & apply_dirichlet_bc,
                 const double         & t = 0);

   private:
      void apply_Dirichlet_BC(SparseMatrix<double> &A,
                              Vector<double>       &b,
                              Vector<double>       &x,
                              const double          t);

      const unsigned int linear_solver_option;

      const ConstraintMatrix * const constraints;

      const DoFHandler<dim> * const dof_handler;

      Function<dim> * const dirichlet_value_function;
};

#include "LinearSolver.cc"
#endif
