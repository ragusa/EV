/**
 * \file NonlinearSolver.h
 * \brief Provides the header for the NonlinearSolver class.
 */

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
 * \brief Class for nonlinear solver.
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

  void initialize(Vector<double> & solution_guess);

  bool update(const SparseMatrix<double> & A, const Vector<double> & b);

protected:

  /** \brief exception for reaching the maximum iteration */
  DeclException1(ExcMaxIterationReached, unsigned int,
    << "Max iteration reached: " << arg1);

  bool check_convergence(const Vector<double> & residual);

  /** \brief Relaxation factor for solution update: \f$\alpha\f$ in
             \f$\mathbf{U}^{(\ell+1)} = \mathbf{U}^{(\ell)}
             + \alpha\Delta\mathbf{U}^{(\ell+1)}\f$ */
  const double relaxation_factor;

  /** \brief Tolerance for convergence of nonlinear system, \f$\epsilon\f$ */
  const double nonlinear_tolerance;

  /** \brief Maximum number of iterations */
  const double iteration_max;

  /** \brief Constraint matrix */
  const ConstraintMatrix & constraints;

  /** \brief Degree of freedom handler, used in linear solver */
  const DoFHandler<dim> & dof_handler;

  /** \brief function for Dirichlet BC values, used in linear solver */
  Function<dim> & dirichlet_value_function;

  /** \brief linear solver */
  LinearSolver<dim> linear_solver;

  /** \brief number of degrees of freedom */
  unsigned int n_dofs;

  /** \brief pointer to solution iterate \f$\mathbf{U}^{(\ell)}\f$ */
  Vector<double> * solution;

  /** \brief change in solution iterate \f$\Delta\mathbf{U}^{(\ell+1)}\f$ */
  Vector<double> solution_change;

  /** \brief The linear residual vector
             \f$\mathbf{r} \equiv \mathbf{b} - \mathbf{A}\mathbf{U}\f$ */
  Vector<double> residual;

  /** \brief Current iteration number */
  unsigned int iteration_number;
};

#include "NonlinearSolver.cc"
#endif
