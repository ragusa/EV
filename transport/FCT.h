#ifndef FCT_cc
#define FCT_cc

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "LinearSolver.h"

using namespace dealii;

/** \brief FCT Class.
 */
template<int dim>
class FCT {
   public:
      FCT(
         const DoFHandler<dim>      &dof_handler,
         const SparseMatrix<double> &lumped_mass_matrix,
         const SparseMatrix<double> &consistent_mass_matrix,
         const LinearSolver<dim>    &linear_solver,
         const SparsityPattern      &sparsity_pattern,
         const std::vector<unsigned int> &dirichlet_nodes,
         const unsigned int         &n_dofs,
         const unsigned int         &dofs_per_cell,
         const bool                 &do_not_limit);
      ~FCT();

      void solve_FCT_system(Vector<double>             &new_solution,
                            const Vector<double>       &old_solution,
                            const SparseMatrix<double> &low_order_ss_matrix,
                            const Vector<double>       &ss_rhs,
                            const double               &dt,
                            const SparseMatrix<double> &low_order_diffusion_matrix,
                            const SparseMatrix<double> &high_order_diffusion_matrix);
      bool check_DMP_satisfied();

   private:
      void compute_bounds(const Vector<double>       &old_solution,
                          const SparseMatrix<double> &low_order_ss_matrix,
                          const Vector<double>       &ss_rhs,
                          const double               &dt);
      void compute_steady_state_bounds(
         const Vector<double>       &old_solution,
         const SparseMatrix<double> &low_order_ss_matrix,
         const Vector<double>       &ss_rhs,
         const double               &dt);
      void compute_flux_corrections(const Vector<double>       &high_order_solution,
                                    const Vector<double>       &old_solution,
                                    const double               &dt,
                                    const SparseMatrix<double> &low_order_diffusion_matrix,
                                    const SparseMatrix<double> &high_order_diffusion_matrix);
      void compute_limiting_coefficients(const Vector<double>       &old_solution,
                                         const SparseMatrix<double> &low_order_ss_matrix,
                                         const Vector<double>       &ss_rhs,
                                         const double               &dt);
      void get_matrix_row(const SparseMatrix<double> &matrix,
                          const unsigned int         &i,
                          std::vector<double>        &row_values,
                          std::vector<unsigned int>  &row_indices,
                          unsigned int               &n_col);
      bool check_max_principle(const Vector<double>       &new_solution,
                               const SparseMatrix<double> &low_order_ss_matrix,
                               const double               &dt);
      void debug_max_principle_low_order(const unsigned int         &i,
                                         const SparseMatrix<double> &low_order_ss_matrix,
                                         const double               &dt);

      const DoFHandler<dim> *dof_handler;

      const SparseMatrix<double> *lumped_mass_matrix;
      const SparseMatrix<double> *consistent_mass_matrix;
      SparsityPattern sparsity_pattern;
      SparseMatrix<double> flux_correction_matrix;
      Vector<double> min_bound;
      Vector<double> max_bound;
      Vector<double> solution_min;
      Vector<double> solution_max;

      SparseMatrix<double> system_matrix;
      Vector<double>       system_rhs;
      Vector<double>       tmp_vector;
      Vector<double>       flux_correction_vector;
      Vector<double>       Q_minus;
      Vector<double>       Q_plus;
      Vector<double>       R_minus;
      Vector<double>       R_plus;

      LinearSolver<dim> linear_solver;

      const std::vector<unsigned int> dirichlet_nodes;

      const unsigned int n_dofs;
      const unsigned int dofs_per_cell;

      const bool do_not_limit;

      bool DMP_satisfied_at_all_steps;
};

#include "FCT.cc"
#endif
