/**
 * \file LinearSolver.h
 * \brief Provides the header for the LinearSolver class.
 */
#ifndef LinearSolver_cc
#define LinearSolver_cc

#include <deal.II/base/function.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include "include/parameters/RunParameters.h"

using namespace dealii;

/**
 * \brief Class for solving linear systems, optionally with Dirichlet BC.
 */
template <int dim>
class LinearSolver
{
public:
  using LinearSolverType = typename RunParameters::LinearSolverType;
  using PreconditionerType = typename RunParameters::PreconditionerType;

  LinearSolver(const RunParameters & parameters,
               const ConstraintMatrix & constraints,
               const DoFHandler<dim> & dof_handler,
               std::shared_ptr<Function<dim>> dirichlet_function,
               const unsigned int & n_components);

  void solve(const SparseMatrix<double> & A,
             Vector<double> & x,
             const Vector<double> & b) const;

  void solve_with_dirichlet(SparseMatrix<double> & A,
                            Vector<double> & x,
                            Vector<double> & b,
                            const bool & dirichlet_bc_apply = true,
                            const double & t = 0.0);

private:
  void apply_dirichlet_bc(SparseMatrix<double> & A,
                          Vector<double> & b,
                          Vector<double> & x,
                          const double & t);

  const LinearSolverType linear_solver_type;

  const PreconditionerType preconditioner_type;

  const unsigned int max_iterations;

  const double linear_tolerance;

  const unsigned int preconditioner_relaxation;

  const bool print_linear_residuals;

  const ConstraintMatrix * const constraints;

  const DoFHandler<dim> * const dof_handler;

  std::shared_ptr<Function<dim>> dirichlet_function;

  const bool have_dirichlet_bc;

  const unsigned int n_components;
};

#include "src/solvers/LinearSolver.cc"
#endif
