#ifndef NonlinearSolver_cc
#define NonlinearSolver_cc

#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include "LinearSolver.h"
#include "TransportParameters.h"

using namespace dealii;

/**
 * Class for nonlinear solver.
 */
template<int dim>
class NonlinearSolver
{
public:

  NonlinearSolver(
    const TransportParameters<dim> & parameters,
    const ConstraintMatrix         & constraints,
    const DoFHandler<dim>          & dof_handler,
    Function<dim>                  & dirichlet_value_function);

  void reset(const Vector<double> & solution_guess);
  bool checkConvergence(const Vector<double> & new_solution);
  Vector<double> getSolution() const;

private:

  /** exception for reaching the maximum iteration */
  DeclException1(ExcMaxIterationReached, unsigned int,
    << "Max iteration reached: " << arg1);

  const double nonlinear_tolerance;
  const double iteration_max;

  const ConstraintMatrix & constraints;
  const DoFHandler<dim> & dof_handler;
  Function<dim> & dirichlet_value_function;
  LinearSolver<dim> linear_solver;

  Vector<double> solution;
  unsigned int iteration_number;
};

#include "NonlinearSolver.cc"
#endif
