#ifndef SSPRKTimeIntegrator_cc
#define SSPRKTimeIntegrator_cc

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include "LinearSolver.h"
#include "TransportParameters.h"

using namespace dealii;

/**
 * Class for implementing SSP Runge-Kutta time integration.
 */
template <int dim>
class SSPRKTimeIntegrator
{
public:
  /** \brief Alias for SSPRK method */
  using SSPRKMethod = typename TransportParameters<dim>::SSPRKMethod;

  SSPRKTimeIntegrator(const SSPRKMethod & ssprk_method,
                      const unsigned int & system_size,
                      const LinearSolver<dim> & linear_solver,
                      const SparsityPattern & sparsity_pattern);

  void initialize_time_step(const Vector<double> & old_solution,
                            const double & time_step_size);

  void step(const SparseMatrix<double> & mass_matrix,
            const SparseMatrix<double> & ss_matrix,
            const Vector<double> & ss_rhs,
            const bool & call_complete_stage_solution);

  void complete_stage_solution();

  void get_stage_solution(const unsigned int & i,
                          Vector<double> & new_solution) const;

  void get_intermediate_solution(Vector<double> & solution) const;

  void set_intermediate_solution(Vector<double> & solution);

  void get_new_solution(Vector<double> & new_solution) const;

  double get_stage_time() const;

  void set_stage_solution(const unsigned int & i,
                          const Vector<double> & solution);

  /** \brief Number of stages in SSPRK method */
  unsigned int n_stages;

private:
  Vector<double> a;
  Vector<double> b;
  Vector<double> c;

  unsigned int n;
  std::vector<Vector<double>> u_stage;

  unsigned int current_stage;

  double t_old;
  double dt;

  LinearSolver<dim> linear_solver;

  SparseMatrix<double> system_matrix;
  Vector<double> system_rhs;
  Vector<double> intermediate_solution;
};

#include "SSPRKTimeIntegrator.cc"
#endif
