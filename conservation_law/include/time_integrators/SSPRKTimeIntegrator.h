/**
 * \file SSPRKTimeIntegrator.h
 * \brief Provides the header for the SSPRKTimeIntegrator class.
 */
#ifndef SSPRKTimeIntegrator_cc
#define SSPRKTimeIntegrator_cc

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include "include/solvers/LinearSolver.h"
#include "include/parameters/RunParameters.h"

using namespace dealii;

/**
 * \brief Class for implementing SSP Runge-Kutta time integration.
 */
template <int dim>
class SSPRKTimeIntegrator
{
public:
  SSPRKTimeIntegrator(const typename RunParameters<
                        dim>::TemporalDiscretization & time_discretization,
                      const unsigned int & system_size,
                      const LinearSolver<dim> & linear_solver,
                      const SparsityPattern & sparsity_pattern);

  void initialize_time_step(const Vector<double> & old_solution,
                            const double & time_step_size);
  void step(const SparseMatrix<double> & mass_matrix,
            const Vector<double> & ss_flux,
            const Vector<double> & ss_rhs,
            const SparseMatrix<double> & diffusion_matrix,
            const bool & call_complete_stage_solution);
  void complete_stage_solution();
  void get_stage_solution(const unsigned int & i,
                          Vector<double> & solution) const;
  void get_intermediate_solution(Vector<double> & solution) const;
  void set_intermediate_solution(const Vector<double> & solution);
  void get_new_solution(Vector<double> & new_solution) const;
  double get_stage_time() const;
  void set_stage_solution(const unsigned int & i,
                          const Vector<double> & solution);

  /** \brief Number of stages */
  unsigned int n_stages;

private:
  /** \brief Stage combination coefficients \f$\alpha_i\f$ */
  std::vector<double> alpha;
  /** \brief Stage combination coefficients \f$\beta_i\f$ */
  std::vector<double> beta;
  /** \brief Stage time coefficients \f$c_i\f$ */
  std::vector<double> c;

  /** \brief Size of linear system */
  unsigned int n;

  /** \brief Vector of the stage solution vectors \f$\hat{\mathrm{U}}^i\f$ */
  std::vector<Vector<double>> u_stage;

  /** \brief Index of current stage */
  unsigned int current_stage;

  /** \brief Old time \f$t^n\f$ */
  double t_old;
  /** \brief Time step size \f$\Delta t\f$ */
  double dt;

  /** \brief Linear solver */
  LinearSolver<dim> linear_solver;

  /** \brief System matrix */
  SparseMatrix<double> system_matrix;
  /** \brief System right hand side vector */
  Vector<double> system_rhs;
  /** \brief Intermediate solution vector \f$\tilde{\mathrm{U}}^i\f$ */
  Vector<double> intermediate_solution;
};

#include "src/time_integrators/SSPRKTimeIntegrator.cc"
#endif
