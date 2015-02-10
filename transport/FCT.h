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
      FCT();
      ~FCT();
      void solve_FCT_system();

   private:
      void compute_bounds();
      void compute_steady_state_bounds();
      void check_bounds();

      SparseMatrix<double> flux_corrections;
      Vector<double> min_bound;
      Vector<double> max_bound;
      Vector<double> solution_min;
      Vector<double> solution_max;

      Vector<double> tmp_vector;
      Vector<double> system_rhs;
};

#include "FCT.cc"
#endif
