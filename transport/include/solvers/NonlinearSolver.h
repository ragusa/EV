/**
 * \file NonlinearSolver.h
 * \brief Provides the header for the NonlinearSolver class.
 */

#ifndef NonlinearSolver_cc
#define NonlinearSolver_cc

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>

#include "include/other/Exceptions.h"
#include "include/solvers/LinearSolver.h"
#include "include/parameters/RunParameters.h"

using namespace dealii;

/**
 * \brief Class for nonlinear solver.
 */
template <int dim>
class NonlinearSolver
{
public:
  NonlinearSolver(const RunParameters & parameters,
                  const LinearSolver<dim> & linear_solver,
                  const ConstraintMatrix & constraints);

  void reinit(const unsigned int & n_dofs);

  void initialize(Vector<double> & solution_guess);

  bool update(const SparseMatrix<double> & A, const Vector<double> & b);

protected:
  bool check_convergence(const Vector<double> & residual);

  /** \brief Relaxation factor for solution update: \f$\alpha\f$ in
             \f$\mathbf{U}^{(\ell+1)} = \mathbf{U}^{(\ell)}
             + \alpha\Delta\mathbf{U}^{(\ell+1)}\f$ */
  const double relaxation_factor;

  /** \brief Tolerance for convergence of nonlinear system, \f$\epsilon\f$ */
  const double nonlinear_tolerance;

  /** \brief Maximum number of iterations */
  const double iteration_max;

  /** \brief pointer to linear solver */
  const LinearSolver<dim> * const linear_solver;

  /** \brief Constraint matrix */
  const ConstraintMatrix & constraints;

  /** \brief Current iteration number */
  unsigned int iteration_number;

  /** \brief number of degrees of freedom */
  unsigned int n_dofs;

  /** \brief pointer to solution iterate \f$\mathbf{U}^{(\ell)}\f$ */
  Vector<double> * solution;

  /** \brief change in solution iterate \f$\Delta\mathbf{U}^{(\ell+1)}\f$ */
  Vector<double> solution_change;

  /** \brief The linear residual vector
             \f$\mathbf{r} \equiv \mathbf{b} - \mathbf{A}\mathbf{U}\f$ */
  Vector<double> residual;

  /** \brief Conditional output stream 1 */
  ConditionalOStream cout1;
};

#include "NonlinearSolver.cc"
#endif
