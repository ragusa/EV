#ifndef FCT_cc
#define FCT_cc

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria.h>
#include "LinearSolver.h"
#include "NonlinearSolver.h"
#include "PostProcessor.h"

using namespace dealii;

/**
 * Class for performing FCT.
 */
template <int dim>
class FCT
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  FCT(const DoFHandler<dim> & dof_handler,
      Triangulation<dim> & triangulation,
      const SparseMatrix<double> & lumped_mass_matrix,
      const SparseMatrix<double> & consistent_mass_matrix,
      const LinearSolver<dim> & linear_solver,
      const SparsityPattern & sparsity_pattern,
      const std::vector<unsigned int> & dirichlet_nodes,
      const unsigned int & n_dofs,
      const unsigned int & dofs_per_cell,
      const bool & do_not_limit,
      const bool & include_analytic_bounds,
      const FESystem<dim> & fe,
      const QGauss<dim> & cell_quadrature,
      const FunctionParser<dim> & cross_section_function,
      FunctionParser<dim> & source_function,
      const bool & source_is_time_dependent,
      const double & theta);

  FCT(const DoFHandler<dim> & dof_handler,
      Triangulation<dim> & triangulation,
      const LinearSolver<dim> & linear_solver,
      const SparsityPattern & sparsity_pattern,
      const std::vector<unsigned int> & dirichlet_nodes,
      const unsigned int & n_dofs,
      const unsigned int & dofs_per_cell,
      const bool & do_not_limit,
      const bool & include_analytic_bounds,
      const FESystem<dim> & fe,
      const QGauss<dim> & cell_quadrature,
      const FunctionParser<dim> & cross_section_function,
      FunctionParser<dim> & source_function);

  void solve_FCT_system_fe(
    Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs,
    const double & dt,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix,
    const double & t_old);

  void compute_bounds_fe(const Vector<double> & old_solution,
                         const SparseMatrix<double> & low_order_ss_matrix,
                         const Vector<double> & ss_rhs,
                         const double & dt,
                         const double & t_old);

  void compute_bounds_theta(const Vector<double> & new_solution,
                            const Vector<double> & old_solution,
                            const SparseMatrix<double> & low_order_ss_matrix,
                            const Vector<double> & ss_rhs_new,
                            const Vector<double> & ss_rhs_old,
                            const double & dt,
                            const double & t_old);

  void compute_bounds_ss(const Vector<double> & solution,
                         const SparseMatrix<double> & low_order_ss_matrix,
                         const Vector<double> & ss_rhs);

  void compute_flux_corrections_fe(
    const Vector<double> & high_order_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix);

  void compute_flux_corrections_theta(
    const Vector<double> & high_order_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix_new,
    const SparseMatrix<double> & high_order_diffusion_matrix_old);

  void compute_flux_corrections_ss(
    const Vector<double> & high_order_solution,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix);

  void compute_limited_flux_bounds_fe(
    const Vector<double> & old_solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs,
    const double & dt);

  void compute_limited_flux_bounds_theta(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs_new,
    const Vector<double> & ss_rhs_old,
    const Vector<double> & cumulative_antidiffusion,
    const double & dt);

  void compute_limited_flux_bounds_ss(
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs);

  void compute_limited_fluxes();

  void compute_limited_flux_matrix();

  bool check_DMP_satisfied();

  bool check_fct_bounds(const Vector<double> & solution) const;

  void output_bounds(PostProcessor<dim> & postprocessor) const;

  Vector<double> get_limited_flux_vector() const;

  void compute_row_sums(const SparseMatrix<double> & matrix,
                        Vector<double> & row_sums) const;

  void subtract_limited_flux_correction_matrix();

protected:
  void get_matrix_row(const SparseMatrix<double> & matrix,
                      const unsigned int & i,
                      std::vector<double> & row_values,
                      std::vector<unsigned int> & row_indices,
                      unsigned int & n_col);

  void debug_max_principle_low_order(
    const unsigned int & i,
    const SparseMatrix<double> & low_order_ss_matrix,
    const double & dt);

  void compute_function_bounds(FunctionParser<dim> & function,
                               const double & time,
                               Vector<double> & min_values,
                               Vector<double> & max_values) const;

  void compute_function_bounds(const FunctionParser<dim> & function,
                               Vector<double> & min_values,
                               Vector<double> & max_values) const;

  void compute_min_max_solution(const Vector<double> & solution,
                                Vector<double> & min_solution,
                                Vector<double> & max_solution) const;

  void compute_and_apply_analytic_bounds(const Vector<double> & solution,
                                         const double & dt,
                                         const double & t_old);

  const DoFHandler<dim> * dof_handler;
  Triangulation<dim> * const triangulation;

  const SparseMatrix<double> * lumped_mass_matrix;
  const SparseMatrix<double> * consistent_mass_matrix;
  SparsityPattern sparsity_pattern;

  /** \brief Matrix of antidiffusive correction fluxes \f$\mathbf{P}\f$ */
  SparseMatrix<double> flux_correction_matrix;

  /** \brief Matrix of limited antidiffusive fluxes
   * \f$\mathbf{L}\odot\mathbf{P}\f$ */
  SparseMatrix<double> limited_flux_matrix;

  Vector<double> solution_min;
  Vector<double> solution_max;
  Vector<double> solution_min_new;
  Vector<double> solution_max_new;
  Vector<double> reaction_min;
  Vector<double> reaction_max;
  Vector<double> source_min;
  Vector<double> source_max;
  Vector<double> lower_bound_analytic;
  Vector<double> upper_bound_analytic;

  SparseMatrix<double> system_matrix;
  Vector<double> system_rhs;
  Vector<double> tmp_vector;
  Vector<double> flux_correction_vector;
  Vector<double> Q_minus;
  Vector<double> Q_plus;
  Vector<double> R_minus;
  Vector<double> R_plus;

  LinearSolver<dim> linear_solver;

  const std::vector<unsigned int> dirichlet_nodes;

  const unsigned int n_dofs;
  const unsigned int dofs_per_cell;

  const bool do_not_limit;

  const bool include_analytic_bounds;

  const FESystem<dim> * const fe;

  const QGauss<dim> * const cell_quadrature;

  const unsigned int n_q_points_cell;

  const FunctionParser<dim> * const cross_section_function;

  FunctionParser<dim> * const source_function;

  const bool source_is_time_dependent;

  bool DMP_satisfied;

  /** \brief Theta parameter for theta time discretization schemes */
  const double theta;
};

#include "FCT.cc"
#endif
