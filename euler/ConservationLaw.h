/**
 * \file ConservationLaw.h
 * \brief Provides the header for the ConservationLaw class.
 */

#ifndef ConservationLaw_h
#define ConservationLaw_h

#include <iostream>
#include <fstream>
#include <algorithm>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include "BoundaryConditions.h"
#include "ConservationLawParameters.h"
#include "DirichletBoundaryConditions.h"
#include "Exceptions.h"
#include "PostProcessor.h"

using namespace dealii;

/**
 * \brief Class for solving a general conservation law.
 *
 * This class solves a general conservation law system of the form:
 * \f[
 *   \frac{\partial\mathbf{u}}{\partial t}
 *   + \nabla \cdot \mathbf{F}(\mathbf{u}) = \mathbf{0},
 * \f]
 * where \f$\mathbf{u}\f$ is the vector of conservation variables and
 * \f$\mathbf{F}(\mathbf{u})\f$ is the flux tensor.
 */
template <int dim>
class ConservationLaw
{
public:
  ConservationLaw(const ConservationLawParameters<dim> & params);
  ~ConservationLaw();
  void run();

protected:
  /** \brief Typedef for cell iterators */
  typedef typename DoFHandler<dim>::active_cell_iterator cell_iterator;

  /** \brief Typedef for cell iterator map to double */
  typedef std::map<cell_iterator, double> cell_map;

  void initialize_system();
  void initialize_runge_kutta();
  void setup_system();
  void compute_error_for_refinement();
  void refine_mesh();
  void update_cell_sizes();
  void assemble_mass_matrix();
  void solve_runge_kutta(PostProcessor<dim> & postprocessor);
  double compute_dt_from_cfl_condition();
  double compute_cfl_number(const double & dt) const;
  void linear_solve(const SparseMatrix<double> & A,
                    const Vector<double> & b,
                    Vector<double> & x);
  void apply_Dirichlet_BC(const double & time);
  void update_viscosities(const double & dt, const unsigned int & n);
  virtual void update_old_low_order_viscosity(
    const bool & using_low_order_scheme);
  void update_max_principle_viscosity();
  void compute_viscous_bilinear_forms();
  void compute_viscous_fluxes();
  void add_maximum_principle_viscosity_bilinear_form(Vector<double> & solution);
  double compute_entropy_normalization(const Vector<double> & solution) const;
  virtual double compute_max_entropy_residual(const Vector<double> & new_solution,
                                              const Vector<double> & old_solution,
                                              const double & dt,
                                              const cell_iterator & cell) const;
  void get_dirichlet_nodes();
  void check_nan();
  bool check_DMP(const unsigned int & n) const;
  void compute_max_principle_quantities();
  void get_matrix_row(const SparseMatrix<double> & matrix,
                      const unsigned int & i,
                      std::vector<double> & row_values,
                      std::vector<unsigned int> & row_indices,
                      unsigned int & n_col);

  virtual void update_entropy_viscosities(const double & dt);
  virtual double compute_max_entropy_jump(const Vector<double> & solution,
                                          const cell_iterator & cell) const;

  /**
   * \brief Computes the lumped mass matrix.
   */
  virtual void assemble_lumped_mass_matrix() = 0;

  /**
   * \brief Computes the flux speed \f$\lambda\f$ at each quadrature point
   *        in domain and finds the max in each cell and the max in the
   *        entire domain.
   */
  virtual void update_flux_speeds() = 0;

  /** \brief Computes the steady-state residual.
   *
   *  This function computes the steady-state residual \f$\mathbf{r}\f$
   *  in the following discretized system:
   *  \f[
   *    \mathbf{M}\frac{d\mathbf{U}}{dt} = \mathbf{r} .
   *  \f]
   *  In general, this steady-state residual consists of an inviscid component
   *  and a viscous component:
   *  \f[
   *    \mathbf{r} = \mathbf{r}^{inviscid} + \mathbf{r}^{viscous} ,
   *  \f]
   *  which are computed as follows:
   *  \f[
   *    r^{inviscid}_i = -\left(\varphi_i,
   *      \nabla\cdot\mathbf{f}(u_h)\right)_\Omega ,
   *  \f]
   *  \f[
   *    r^{viscous}_i = -\sum\limits_{K\subset S_i}\nu_K\sum\limits_j
   *      U_j b_K(\varphi_i, \varphi_j) .
   *  \f]
   *
   *  \param[out] r steady-state residual \f$\mathbf{r}\f$
   */
  virtual void compute_ss_residual(Vector<double> & r) = 0;

  // virtual void compute_ss_jacobian() = 0;

  // void compute_tr_residual(unsigned int i, double dt);

  /**
   * \brief Computes entropy for each quadrature point on a cell or face.
   *
   * \param[in] solution solution
   * \param[in] fe_values FE values, either for a cell or a face
   * \param[out] entropy entropy values at each quadrature point on
   *             cell or face
   */
  virtual void compute_entropy(const Vector<double> & solution,
                               const FEValuesBase<dim> & fe_values,
                               Vector<double> & entropy) const = 0;

  /**
   * \brief Computes divergence of entropy flux at each quadrature point in cell.
   *
   * \param[in] solution solution
   * \param[in] fe_values FEValues object
   * \param[out] divergence divergence of entropy flux at each quadrature
   *             point on cell
   */
  virtual void compute_divergence_entropy_flux(const Vector<double> &,
                                               const FEValuesBase<dim> &,
                                               Vector<double> &) const
  {
  }

  /**
   * \brief Outputs additional quantities other than the solution variables.
   *
   * Derived classes define this function if any additional quantities are
   * desired to be output.
   *
   * \param[in] postprocessor the post-processor object
   */
  virtual void output_results(PostProcessor<dim> & postprocessor) const;

  /**
   * \brief Returns the names of each component.
   *
   * \return vector of names of each component
   */
  virtual std::vector<std::string> get_component_names() = 0;

  /**
   * \brief Returns the interpretations for each component.
   *
   * This function returns the interpretation of each component,
   * i.e., whether each component is a scalar or a component of
   * a vector.
   *
   * \return data component interpretations
   */
  virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_component_interpretations() = 0;

  /**
   * \brief Creates the domain, computes volume, defines initial
   *        conditions, and defines boundary conditions and exact
   *        solution if it exists.
   */
  virtual void define_problem() = 0;

  /**
   * \brief Performs additional setup required by a derived class.
   */
  virtual void perform_additional_setup() {}
  /** \brief name of problem */
  std::string problem_name;

  /** \brief input parameters for conservation law */
  const ConservationLawParameters<dim> parameters;

  /** \brief number of components in the system */
  const unsigned int n_components;

  /** \brief triangulation */
  Triangulation<dim> triangulation;
  /** \brief mapping */
  const MappingQ1<dim> mapping;

  /** \brief finite element system */
  const FESystem<dim> fe;
  /** \brief DoFs per cell */
  const unsigned int dofs_per_cell;
  /** \brief faces per cell */
  const unsigned int faces_per_cell;
  /** \brief Number of degrees of freedom */
  unsigned int n_dofs;
  /** \brief Degree of freedom handler */
  DoFHandler<dim> dof_handler;
  /** \brief constraint matrix */
  ConstraintMatrix constraints;

  /** \brief number of quadrature points in each dimension */
  const unsigned int n_q_points_per_dim;
  /** \brief quadrature formula for cells */
  const QGauss<dim> cell_quadrature;
  /** \brief quadrature formula for faces */
  const QGauss<dim - 1> face_quadrature;
  /** \brief number of quadrature points per cell */
  const unsigned int n_q_points_cell;
  /** \brief number of quadrature points per face */
  const unsigned int n_q_points_face;

  /** \brief solution of current time step */
  Vector<double> new_solution;
  /** \brief solution of previous time step */
  Vector<double> old_solution;

  /** \brief exact solution of current time step */
  Vector<double> exact_solution;

  /** \brief solution step in Newton loop; vector for temporary storage */
  Vector<double> solution_step;
  /** \brief current time */
  double current_time;
  /** \brief old time (beginning of each time step) */
  double old_time;
  /** \brief system right-hand side */
  Vector<double> system_rhs;

  /** \brief constrained sparsity pattern */
  SparsityPattern constrained_sparsity_pattern;
  /** \brief unconstrained sparsity pattern */
  SparsityPattern unconstrained_sparsity_pattern;
  /** \brief consistent mass matrix */
  SparseMatrix<double> consistent_mass_matrix;
  /** \brief lumped mass matrix */
  SparseMatrix<double> lumped_mass_matrix;
  /** \brief system matrix; Jacobian matrix */
  SparseMatrix<double> system_matrix;
  /**
   * \brief Matrix for the viscous bilinear forms, to be used in computing
   *        viscosity.
   *
   * Each element of this matrix, \f$B_{i,j}\f$ is the following:
   * \f[
   *   B_{i,j} = \sum_{K\subset S_{i,j}}b_K(\varphi_i,\varphi_j)
   * \f]
   */
  SparseMatrix<double> viscous_bilinear_forms;

  /** \brief Viscous fluxes */
  SparseMatrix<double> viscous_fluxes;

  /** \brief number of boundaries */
  unsigned int n_boundaries;
  /** \brief type of boundary conditions */
  std::string boundary_conditions_type;
  /** \brief boundary conditions */
  std::shared_ptr<BoundaryConditions<dim>> boundary_conditions;
  /** \brief option to use exact solution function as Dirichlet boundary
   *         conditions */
  bool use_exact_solution_as_dirichlet_bc;
  /** \brief vector of Dirichlet BC function strings, which will be parsed */
  std::vector<std::vector<std::string>> dirichlet_function_strings;
  /** \brief vector of Dirichlet BC functions created from parsed strings */
  std::vector<FunctionParser<dim> *> dirichlet_function;

  /** \brief initial conditions function strings for each component, which will
   *         be parsed */
  std::vector<std::string> initial_conditions_strings;
  /** \brief initial conditions functions */
  FunctionParser<dim> initial_conditions_function;

  /** \brief constants for function parsers */
  std::map<std::string, double> constants;

  /** \brief option if the problem has an exact solution provided */
  bool has_exact_solution;
  /** \brief Exact solution function strings for each component */
  std::vector<std::string> exact_solution_strings;
  /** \brief Exact solution function */
  std::shared_ptr<Function<dim>> exact_solution_function;
  /** \brief Convergence table for errors computed in each refinement cycle **/
  ConvergenceTable convergence_table;

  /** \brief Vector of component names */
  std::vector<std::string> component_names;
  /** \brief Vector of data component interpretations (scalar or vector) */
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretations;

  /** \brief Minimum cell diameter */
  double minimum_cell_diameter;
  /** \brief Maximum flux speed in domain */
  double max_flux_speed;
  /** \brief Domain volume */
  double domain_volume;

  // maps
  cell_map cell_diameter;
  cell_map max_flux_speed_cell;
  cell_map viscosity;
  cell_map first_order_viscosity;
  cell_map entropy_viscosity;

  /** \brief Flag to signal being in last adaptive refinement cycle */
  bool in_final_cycle;
  /** \brief Estimation of error per cell for adaptive mesh refinement */
  Vector<float> estimated_error_per_cell;

  /**
   * \brief Structure for Runge-Kutta constants and memory for stage
   *        steady-state residuals.
   */
  struct RungeKuttaParameters
  {
    int s;
    std::vector<Vector<double>> a;
    std::vector<double> b;
    std::vector<double> c;

    bool solution_computed_in_last_stage;
    bool is_explicit;

    std::vector<Vector<double>> f;
  };

  /** \brief Runge-Kutta parameters */
  RungeKuttaParameters rk;

  /** \brief Minimum values for use in DMP */
  Vector<double> min_values;
  /** \brief Maximum values for use in DMP */
  Vector<double> max_values;

  /** \brief Vector of degrees of freedom subject to Dirichlet boundary
   *         conditions. */
  std::vector<unsigned int> dirichlet_nodes;

  /** \brief Default end time for test problem */
  double default_end_time;
  /** \brief Flag to signal that test problem has a default end time */
  bool has_default_end_time;
  /** \brief Chosen end time */
  double end_time;
};

#include "ConservationLaw.cc"

#endif
