/**
 * \file FCT.h
 * \brief Provides the header for the FCT class.
 */

#ifndef FCT_h
#define FCT_h

#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_accessor.h>
#include "include/parameters/ConservationLawParameters.h"
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
  using AntidiffusionType =
    typename ConservationLawParameters<dim>::AntidiffusionType;

  using FCTSynchronizationType =
    typename ConservationLawParameters<dim>::FCTSynchronizationType;

  using FCTVariablesType =
    typename ConservationLawParameters<dim>::FCTVariablesType;

  FCT(const ConservationLawParameters<dim> & parameters,
      const DoFHandler<dim> & dof_handler,
      const Triangulation<dim> & triangulation,
      const SparseMatrix<double> & lumped_mass_matrix,
      const SparseMatrix<double> & consistent_mass_matrix,
      const std::shared_ptr<StarState<dim>> & star_state,
      const LinearSolver<dim> & linear_solver,
      const SparsityPattern & sparsity_pattern,
      const std::vector<unsigned int> & dirichlet_nodes,
      const unsigned int & n_components,
      const unsigned int & dofs_per_cell);

  void solve_fct_system(Vector<double> & new_solution,
                        const Vector<double> & old_solution,
                        const Vector<double> & ss_flux,
                        const Vector<double> & ss_reaction,
                        const Vector<double> & ss_rhs,
                        const double & dt,
                        const SparseMatrix<double> & low_order_diffusion_matrix,
                        const SparseMatrix<double> & high_order_diffusion_matrix);

  bool check_DMP_satisfied();

  /*
    void output_bounds(const PostProcessor<dim> & postprocessor,
                       const std::string & description_string) const;
  */

  virtual void compute_bounds(const Vector<double> & old_solution,
                              const Vector<double> & ss_reaction,
                              const Vector<double> & ss_rhs,
                              const double & dt);

  void output_limiter_matrix() const;

private:
  void compute_flux_corrections(
    const Vector<double> & high_order_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix);

  void compute_limiting_coefficients_zalesak(
    const Vector<double> & old_solution,
    const Vector<double> & ss_flux,
    const Vector<double> & ss_rhs,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const double & dt);

  void compute_limited_flux_correction_vector();

  void compute_full_flux_correction();

  void get_matrix_row(const SparseMatrix<double> & matrix,
                      const unsigned int & i,
                      std::vector<double> & row_values,
                      std::vector<unsigned int> & row_indices,
                      unsigned int & n_col);

  void synchronize_min();

  void check_limited_flux_correction_sum(
    const Vector<double> & flux_correction_vector);

  /*
  bool check_max_principle(const Vector<double> & new_solution,
                           const SparseMatrix<double> & low_order_ss_matrix,
                           const double & dt);
  void debug_max_principle_low_order(
    const unsigned int & i,
    const SparseMatrix<double> & low_order_ss_matrix,
    const double & dt);
*/
  const FE_Q<dim> fe_scalar;

  /** \brief scalar sparsity pattern */
  SparsityPattern scalar_sparsity_pattern;

  const DoFHandler<dim> * dof_handler;

  DoFHandler<dim> dof_handler_scalar;

  const SparseMatrix<double> * const lumped_mass_matrix;
  const SparseMatrix<double> * const consistent_mass_matrix;
  SparsityPattern sparsity_pattern;
  SparseMatrix<double> limiter_matrix;
  SparseMatrix<double> flux_correction_matrix;
  Vector<double> solution_min;
  Vector<double> solution_max;

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

  const AntidiffusionType antidiffusion_type;

  const FCTSynchronizationType synchronization_type;

  const FCTVariablesType fct_variables_type;

  const bool use_star_states_in_fct_bounds;

  bool DMP_satisfied_at_all_steps;
};

#include "src/fct/FCT.cc"
#endif
