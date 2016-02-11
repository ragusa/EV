/**
 * \file ShallowWaterFCT.h
 * \brief Provides the header for the ShallowWaterFCT class.
 */

#ifndef ShallowWaterFCT_h
#define ShallowWaterFCT_h

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_accessor.h>
#include "include/parameters/RunParameters.h"
#include "include/postprocessing/PostProcessor.h"
#include "include/riemann/StarState.h"
#include "include/solvers/LinearSolver.h"

using namespace dealii;

/**
 * \brief Class for performing FCT for the shallow water equations.
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
class ShallowWaterFCT : public FCT<dim>
{
public:
  using AntidiffusionType =
    typename RunParameters<dim>::AntidiffusionType;

  using FCTSynchronizationType =
    typename RunParameters<dim>::FCTSynchronizationType;

  ShallowWaterFCT(const RunParameters<dim> & parameters,
                  const DoFHandler<dim> & dof_handler,
                  const SparseMatrix<double> & lumped_mass_matrix,
                  const SparseMatrix<double> & consistent_mass_matrix,
                  const std::shared_ptr<StarState<dim>> & star_state,
                  const LinearSolver<dim> & linear_solver,
                  const SparsityPattern & sparsity_pattern,
                  const std::vector<unsigned int> & dirichlet_nodes,
                  const unsigned int & n_components,
                  const unsigned int & dofs_per_cell);

  void compute_bounds(const Vector<double> & old_solution,
                      const Vector<double> & ss_reaction,
                      const Vector<double> & ss_rhs,
                      const double & dt) override;

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
};

#include "src/fct/ShallowWaterFCT.cc"
#endif
