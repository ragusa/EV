/**
 * \file FCT.h
 * \brief Provides the header for the FCT class.
 */

#ifndef FCT_h
#define FCT_h

#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_accessor.h>
#include "include/parameters/RunParameters.h"
#include "include/postprocessing/PostProcessor.h"
#include "include/riemann/StarState.h"
#include "include/solvers/LinearSolver.h"

using namespace dealii;

/**
 * \brief Class for performing FCT.
 *
 * This class allows the FCT system to be solved for the case of Explicit Euler:
 * \f[
 *   \mathrm{M}^L\mathrm{U}^{n+1} = \mathrm{M}^L\mathrm{U}^n + \Delta t\left(
 *     \mathrm{b}^n - \mathbf{C}\mathbf{F}(\mathrm{U}^n)
 *     - \mathrm{D}^L(\mathrm{U}^n)\mathrm{U}^n - \bar{\mathrm{p}}\right) \,,
 * \f]
 * where \f$\bar{\mathrm{p}}\f$ is the limited flux sum vector. The limiting
 * coefficients enforce the following bounds:
 * \f[
 *   W_i^-(\mathrm{U}^n) \leq U_i^{n+1} \leq W_i^+(\mathrm{U}^n) \,,
 * \f]
 * where the bounds are the following:
 * \f[
 *   W_i^-(\mathrm{U}^n) = U_{min,i}^n\left(1
 *     - \frac{\Delta t}{M^L_{i,i}}\sigma_i\right)
 *     + \frac{\Delta t}{M^L_{i,i}}b_i^n \,,
 * \f]
 * \f[
 *   W_i^+(\mathrm{U}^n) = U_{max,i}^n\left(1
 *     - \frac{\Delta t}{M^L_{i,i}}\sigma_i\right)
 *     + \frac{\Delta t}{M^L_{i,i}}b_i^n \,,
 * \f]
 * where \f$\sigma_i\f$ is the reaction vector:
 * \f[
 *   \sigma_i = \int\limits_{S_i}\varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *     = \sum\limits_j\int\limits_{S_{i,j}}
 *     \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \,,
 * \f]
 * where \f$\sigma(\mathbf{x})\f$ is the reaction coefficient of the
 * conservation law equation when it is put in the form
 * \f[
 *   \frac{\partial\mathbf{u}}{\partial t}
 *   + \nabla \cdot \mathbf{F}(\mathbf{u}) + \sigma(\mathbf{x})\mathbf{u}
 *   = \mathbf{0} \,.
 * \f]
 */
template <int dim>
class FCT
{
public:
  using FCTBoundsType = typename RunParameters<dim>::FCTBoundsType;

  using AntidiffusionType = typename RunParameters<dim>::AntidiffusionType;

  using FCTSynchronizationType =
    typename RunParameters<dim>::FCTSynchronizationType;

  using FCTLimitationType = typename RunParameters<dim>::FCTLimitationType;

  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  FCT(const RunParameters<dim> & parameters,
      const DoFHandler<dim> & dof_handler,
      const Triangulation<dim> & triangulation,
      const SparseMatrix<double> & lumped_mass_matrix,
      const SparseMatrix<double> & consistent_mass_matrix,
      const std::shared_ptr<StarState<dim>> & star_state,
      const LinearSolver<dim> & linear_solver,
      const SparsityPattern & sparsity_pattern,
      const std::vector<unsigned int> & dirichlet_nodes,
      const unsigned int & n_components,
      const unsigned int & dofs_per_cell,
      const std::vector<std::string> & component_names,
      const bool & use_star_states_in_fct_bounds);

  void reinitialize(const SparsityPattern & sparsity_pattern);

  void solve_fct_system(Vector<double> & new_solution,
                        const Vector<double> & old_solution,
                        const Vector<double> & ss_flux,
                        const Vector<double> & ss_reaction,
                        const Vector<double> & ss_rhs,
                        const double & dt,
                        const SparseMatrix<double> & low_order_diffusion_matrix,
                        const SparseMatrix<double> & high_order_diffusion_matrix);

  bool check_fct_bounds_satisfied(const Vector<double> & new_solution) const;

  void compute_min_and_max_of_solution(const Vector<double> & solution,
                                       Vector<double> & min_values,
                                       Vector<double> & max_values) const;

  void compute_solution_bounds(const Vector<double> & old_solution,
                               const Vector<double> & ss_reaction,
                               const Vector<double> & ss_rhs,
                               const double & dt);

  void compute_solution_bounds_characteristic(const Vector<double> & old_solution,
                                              const Vector<double> & ss_reaction,
                                              const Vector<double> & ss_rhs,
                                              const double & dt);

  void output_limiter_matrix() const;

  void output_bounds_transient(PostProcessor<dim> & postprocessor,
                               const double & time);

protected:
  void compute_flux_corrections(
    const Vector<double> & high_order_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix);

  void compute_antidiffusion_bounds(
    const Vector<double> & old_solution,
    const Vector<double> & ss_flux,
    const Vector<double> & ss_rhs,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const double & dt);

  void compute_antidiffusion_bounds_characteristic(
    const Vector<double> & old_solution, const double & dt);

  void compute_limiting_coefficients_zalesak();

  void compute_limited_flux_correction_vector();

  void compute_full_flux_correction_vector();

  void get_matrix_row(const SparseMatrix<double> & matrix,
                      const unsigned int & i,
                      std::vector<double> & row_values,
                      std::vector<unsigned int> & row_indices,
                      unsigned int & n_col);

  void synchronize_min();

  virtual FullMatrix<double> compute_transformation_matrix(
    const Vector<double> & solution) const;

  virtual FullMatrix<double> compute_transformation_matrix_inverse(
    const Vector<double> & solution) const;

  void transform_vector(const Vector<double> & solution,
                        const Vector<double> & vector_original,
                        Vector<double> & vector_transformed) const;

  void transform_matrix(const Vector<double> & solution,
                        const SparseMatrix<double> & matrix_original,
                        SparseMatrix<double> & matrix_transformed) const;

  void create_dof_indices_lists();

  void check_conservation(const Vector<double> & flux_correction_vector);

  const FE_Q<dim> fe_scalar;

  /** \brief scalar sparsity pattern */
  SparsityPattern scalar_sparsity_pattern;

  const DoFHandler<dim> * dof_handler;

  DoFHandler<dim> dof_handler_scalar;

  const SparseMatrix<double> * const lumped_mass_matrix;
  const SparseMatrix<double> * const consistent_mass_matrix;
  SparsityPattern sparsity_pattern;
  SparseMatrix<double> limiter_matrix;

  /** \brief matrix of antidiffusive fluxes \f$\mathbf{P}\f$ */
  SparseMatrix<double> flux_correction_matrix;

  /** \brief matrix of transformed antidiffusive fluxes \f$\hat{\mathbf{P}}\f$ */
  SparseMatrix<double> flux_correction_matrix_transformed;

  Vector<double> solution_min;
  Vector<double> solution_max;

  /** \brief old solution vector in characteristic variables \hat{\mathbf{U}}^n */
  Vector<double> old_solution_characteristic;

  SparseMatrix<double> system_matrix;
  Vector<double> system_rhs;
  Vector<double> tmp_vector;
  Vector<double> flux_correction_vector;
  Vector<double> Q_minus;
  Vector<double> Q_plus;
  Vector<double> R_minus;
  Vector<double> R_plus;

  std::shared_ptr<StarState<dim>> star_state;

  LinearSolver<dim> linear_solver;

  const SparsityPattern * const sparsity;

  const std::vector<unsigned int> dirichlet_nodes;

  const unsigned int n_dofs;

  const unsigned int n_components;

  const unsigned int n_dofs_scalar;

  const unsigned int dofs_per_cell;

  const unsigned int dofs_per_cell_per_component;

  const FCTBoundsType fct_bounds_type;

  const AntidiffusionType antidiffusion_type;

  const FCTSynchronizationType synchronization_type;

  const FCTLimitationType fct_limitation_type;

  const bool use_star_states_in_fct_bounds;

  bool fct_bounds_satisfied_at_all_steps;

  unsigned int bounds_transient_file_index;

  std::vector<std::string> lower_bound_component_names;

  std::vector<std::string> upper_bound_component_names;

  /** \brief vector of times and corresponding lower bound file names */
  std::vector<std::pair<double, std::string>> times_and_lower_bound_filenames;

  /** \brief vector of times and corresponding upper bound file names */
  std::vector<std::pair<double, std::string>> times_and_upper_bound_filenames;

  /** \brief list of node index of each degree of freedom */
  std::vector<unsigned int> node_indices;

  /** \brief list of component index of each degree of freedom */
  std::vector<unsigned int> component_indices;

  /** \brief list of degree of freedom indices for each node */
  std::vector<std::vector<unsigned int>> node_dof_indices;

  /** \brief number of nodes */
  unsigned int n_nodes;
};

#include "src/fct/FCT.cc"
#endif
