/**
 * \file FCT.h
 * \brief Provides the header for the FCT class.
 */
#ifndef FCT_h
#define FCT_h

#include <sstream>
#include <string>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "include/fct/MultipassLimiter.h"
#include "include/fct/OnesLimiter.h"
#include "include/fct/ZalesakLimiter.h"
#include "include/fct/ZeroesLimiter.h"
#include "include/parameters/RunParameters.h"
#include "include/postprocessing/PostProcessor.h"

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
 *
 * \note This class and its children must be instantiated in each refinement
 *       cycle because there is currently no ability to reinitialize to a
 *       new system size.
 */
template <int dim>
class FCT
{
public:
  /** \brief alias for limiter option */
  using LimiterOption = typename RunParameters::LimiterOption;

  FCT(const RunParameters & run_parameters,
      const DoFHandler<dim> & dof_handler,
      const FESystem<dim> & fe,
      const std::map<unsigned int, double> & dirichlet_values);

  virtual bool check_bounds(const Vector<double> & new_solution) const = 0;

  void output_bounds(PostProcessor<dim> & postprocessor) const;

  void output_bounds_transient(PostProcessor<dim> & postprocessor,
                               const double & time);

  virtual Vector<double> get_lower_solution_bound() const = 0;

  virtual Vector<double> get_upper_solution_bound() const = 0;

  void output_limiter_matrix() const;

protected:
  void compute_row_sum_vector(const SparseMatrix<double> & matrix,
                              Vector<double> & row_sum_vector) const;

  std::vector<std::string> create_filter_string_list(
    std::string filter_sequence_string);

  /** \brief run parameters */
  const RunParameters * const run_parameters;

  /** \brief map of DoF indices to Dirichlet BC values */
  const std::map<unsigned int, double> * const dirichlet_values;

  /** \brief limiter */
  std::shared_ptr<Limiter<dim>> limiter;

  /** \brief sparsity pattern */
  std::shared_ptr<SparsityPattern> sparsity_pattern;

  /** \brief limiter matrix \f$\mathbf{L}\f$ */
  SparseMatrix<double> limiter_matrix;

  /** \brief antidiffusion matrix \f$\mathbf{P}\f$ */
  SparseMatrix<double> antidiffusion_matrix;

  /** \brief limited antidiffusion matrix \f$\bar{\mathbf{P}}\f$ */
  SparseMatrix<double> limited_antidiffusion_matrix;

  /** \brief cumulative antidiffusion vector \f$\bar{\mathbf{p}}\f$ */
  Vector<double> cumulative_antidiffusion;

  /** \brief degree of freedom handler */
  const DoFHandler<dim> * const dof_handler;

  /** \brief finite element system */
  const FESystem<dim> * const fe;

  /** \brief number of degrees of freedom */
  const unsigned int n_dofs;

  /** \brief comma-delimited string of filter identifier strings */
  std::string filter_sequence_string;

  /** \brief number of FCT filters */
  unsigned int n_filters;

  /** \brief index for bounds transient file */
  unsigned int bounds_transient_file_index;

  /** \brief names of lower bounds of each component */
  std::vector<std::string> lower_bound_component_names;

  /** \brief names of upper bounds of each component */
  std::vector<std::string> upper_bound_component_names;

  /** \brief vector of times and corresponding lower bound file names */
  std::vector<std::pair<double, std::string>> times_and_lower_bound_filenames;

  /** \brief vector of times and corresponding upper bound file names */
  std::vector<std::pair<double, std::string>> times_and_upper_bound_filenames;
};

#include "src/fct/FCT.cc"

#endif
