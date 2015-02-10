#ifndef PostProcessor_cc
#define PostProcessor_cc

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "LinearSolver.h"

using namespace dealii;

/** \brief Class for implementing SSP Runge-Kutta time integration
 */
template<int dim>
class PostProcessor {
   public:
      PostProcessor(const SparsistyPattern &sparsity_pattern);
      ~PostProcessor();
      void solve_PostProcessor_system();

   private:
      void compute_bounds();
      void compute_steady_state_bounds();
      void check_bounds();

      SparsityPattern sparsity_pattern;
      SparseMatrix<double> flux_corrections;
      Vector<double> min_bound;
      Vector<double> max_bound;
      Vector<double> solution_min;
      Vector<double> solution_max;

      Vector<double> tmp_vector;
      Vector<double> system_rhs;
};

#include "PostProcessor.cc"
#endif
