#ifndef SSPRungeKuttaTimeIntegrator_cc
#define SSPRungeKuttaTimeIntegrator_cc

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "LinearSolver.h"

using namespace dealii;

/** \brief Class for implementing SSP Runge-Kutta time integration
 */
template<int dim>
class SSPRungeKuttaTimeIntegrator {
   public:
      SSPRungeKuttaTimeIntegrator(const unsigned int      &n_stages,
                                  const unsigned int      &system_size,
                                  const LinearSolver<dim> &linear_solver,
                                  const SparsityPattern   &sparsity_pattern);
      ~SSPRungeKuttaTimeIntegrator();

      void initialize_time_step(const Vector<double> &old_solution,
                                const double         &time_step_size);
      void advance_stage(const SparseMatrix<double> &mass_matrix,
                         const SparseMatrix<double> &ss_matrix,
                         const Vector<double>       &ss_rhs);
      void get_stage_solution(const unsigned int &i, Vector<double> &new_solution) const;
      void get_new_solution(Vector<double> &new_solution) const;
      double get_stage_time() const;
      void set_stage_solution(const unsigned int &i, const Vector<double> &solution);

      unsigned int n_stages;

   private:
      Vector<double> a;
      Vector<double> b;
      Vector<double> c;

      unsigned int n;
      std::vector<Vector<double> > u_stage;

      unsigned int current_stage;

      double t_old;
      double dt;

      LinearSolver<dim> linear_solver;

      SparseMatrix<double> system_matrix;
      Vector<double> system_rhs;
      Vector<double> tmp_vector;
};

#include "SSPRungeKuttaTimeIntegrator.cc"
#endif
