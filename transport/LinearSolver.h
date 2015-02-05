#ifndef LinearSolver_cc
#define LinearSolver_cc

#include <deal.II/base/function.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

/** \brief Class for solving linear systems.
 */
template<int dim>
class LinearSolver {
   public:
      LinearSolver(const unsigned int     &linear_solver_option,
                   const ConstraintMatrix &constraints,
                   const DoFHandler<dim>  &dof_handler,
                   const Function<dim>    &dirichlet_value_function);
      ~LinearSolver();

      void solve(SparseMatrix<double> &A,
                 Vector<double>       &b,
                 Vector<double>       &x) const;

   private:
      void apply_Dirichlet_BC(SparseMatrix<double> &A,
                              Vector<double>       &b,
                              Vector<double>       &x) const;

      const unsigned int linear_solver_option;

      const ConstraintMatrix constraints;

      const DoFHandler<dim> *dof_handler;

      const Function<dim> *dirichlet_value_function;
};

#include "LinearSolver.cc"
#endif
