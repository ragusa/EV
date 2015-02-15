#ifndef FCT_cc
#define FCT_cc

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "LinearSolver.h"

using namespace dealii;

/** \brief Class for implementing SSP Runge-Kutta time integration
 */
template<int dim>
class FCT {
   public:
      FCT(const SparsistyPattern &sparsity_pattern,
          const LinearSolver<dim> &linear_solver);
      ~FCT();
      void solve_FCT_system();

   private:
      void compute_bounds();
      void compute_steady_state_bounds();
      //void check_bounds();

      const SparseMatrix<double> *lumped_mass_matrix;
      const SparseMatrix<double> *consistent_mass_matrix;
      SparsityPattern sparsity_pattern;
      SparseMatrix<double> flux_corrections;
      Vector<double> min_bound;
      Vector<double> max_bound;
      Vector<double> solution_min;
      Vector<double> solution_max;

      SparseMatrix<double> system_matrix;
      Vector<double>       system_rhs;
      Vector<double>       tmp_vector;

      LinearSolver<dim> linear_solver;
};

#include "FCT.cc"
#endif
