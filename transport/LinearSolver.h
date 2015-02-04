#ifndef LinearSolver_cc
#define LinearSolver_cc

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

/** \brief Class for solving linear systems.
 */
class LinearSolver {
   public:
      LinearSolver(const unsigned int     &linear_solver_option,
                   const ConstraintMatrix &constraints);
      ~LinearSolver();

      void solve(const SparseMatrix<double> &A,
                 const Vector<double>       &b,
                 Vector<double>             &x) const;

   private:
      unsigned int linear_solver_option;
      ConstraintMatrix constraints;
};

#include "LinearSolver.cc"
#endif
