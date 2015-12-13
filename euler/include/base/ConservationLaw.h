/**
 * \file ConservationLaw.h
 * \brief Provides the header for the ConservationLaw class.
 */

#ifndef ConservationLaw_h
#define ConservationLaw_h

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/timer.h>
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
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>

#include "include/bc/BoundaryConditions.h"
#include "include/bc/DirichletBoundaryConditions.h"
#include "include/diffusion/ArtificialDiffusion.h"
#include "include/diffusion/NoDiffusion.h"
#include "include/diffusion/GraphTheoreticDiffusion.h"
#include "include/diffusion/LaplacianDiffusion.h"
#include "include/entropy/Entropy.h"
#include "include/other/Exceptions.h"
#include "include/parameters/ConservationLawParameters.h"
#include "include/postprocessing/PostProcessor.h"
#include "include/viscosity/ConstantViscosity.h"
#include "include/viscosity/DomainInvariantViscosity.h"
#include "include/viscosity/EntropyViscosity.h"
#include "include/viscosity/HighOrderViscosity.h"
#include "include/viscosity/LowOrderViscosity.h"
#include "include/viscosity/MaxWaveSpeed.h"
#include "include/viscosity/Viscosity.h"
#include "include/viscosity/ViscosityMultiplier.h"

#ifdef IS_PARALLEL
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/sparsity_tools.h>
#endif

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
#ifdef IS_PARALLEL
  /** \brief Typedef for vector */
  typedef LinearAlgebraPETSc::MPI::Vector LocalVector;

  /** \brief Typedef for sparse matrix */
  typedef LinearAlgebraPETSc::MPI::SparseMatrix LocalMatrix;

  /** \brief Typedef for triangulation */
  typedef parallel::distributed::Triangulation<dim> LocalTriangulation;
#else
  /** \brief Typedef for vector */
  typedef Vector<double> LocalVector;

  /** \brief Typedef for sparse matrix */
  typedef SparseMatrix<double> LocalMatrix;

  /** \brief Typedef for triangulation */
  typedef Triangulation<dim> LocalTriangulation;
#endif

  /** \brief Typedef for cell iterator */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  /** \brief Typedef for face iterator */
  typedef typename DoFHandler<dim>::active_face_iterator Face;

  /** \brief Typedef for cell iterator map to double */
  typedef std::map<Cell, double> CellMap;

  using TimeStepSizeMethod =
    typename ConservationLawParameters<dim>::TimeStepSizeMethod;
  using TemporalIntegrator =
    typename ConservationLawParameters<dim>::TemporalIntegrator;
  using TemporalDiscretization =
    typename ConservationLawParameters<dim>::TemporalDiscretization;
  using ViscosityType = typename ConservationLawParameters<dim>::ViscosityType;
  using DiffusionType = typename ConservationLawParameters<dim>::DiffusionType;
  using LinearSolverType =
    typename ConservationLawParameters<dim>::LinearSolverType;

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
  void check_nan();
  void output_viscosity(PostProcessor<dim> & postprocessor,
                        const bool & is_transient = false,
                        const double & time = 0.0);

  /*
    virtual void update_entropy_viscosities(const double & dt);
    virtual double compute_max_entropy_jump(const Vector<double> & solution,
                                            const Cell & cell) const;
  */

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
   *  \param[in] dt time step size \f$\Delta t\f$
   *  \param[out] r steady-state residual \f$\mathbf{r}\f$
   */
  virtual void compute_ss_residual(const double & dt, Vector<double> & r) = 0;

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
   * \brief Creates an auxiliary post-processor object and returns the pointer.
   *
   * This default version returns a null pointer.
   *
   * \return pointer to created auxiliary post-processor object
   */
  virtual std::shared_ptr<DataPostprocessor<dim>> create_auxiliary_postprocessor()
    const;

  /**
   * \brief Creates a viscosity multiplier object and returns the pointer.
   *
   * \return pointer to created viscosity multiplier object
   */
  virtual std::shared_ptr<ViscosityMultiplier<dim>> create_viscosity_multiplier()
    const;

  /**
   * \brief Creates an entropy object and returns the pointer.
   *
   * \return pointer to created entropy object
   */
  virtual std::shared_ptr<Entropy<dim>> create_entropy() const = 0;

  /**
   * \brief Creates a max wave speed object and returns the pointer.
   *
   * \return pointer to created max wave speed object
   */
  virtual std::shared_ptr<MaxWaveSpeed<dim>> create_max_wave_speed() const = 0;

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
   * \brief Returns vectors of the scalar and vector FE extractors.
   *
   * \param[out] vector of the scalar extractors
   * \param[out] vector of the vector extractors
   */
  virtual void get_fe_extractors(
    std::vector<FEValuesExtractors::Scalar> & scalar_extractors,
    std::vector<FEValuesExtractors::Vector> & vector_extractors) const = 0;

  /**
   * \brief Creates the domain, computes volume, defines initial
   *        conditions, and defines boundary conditions and exact
   *        solution if it exists.
   */
  virtual void define_problem() = 0;

  /**
   * \brief Performs non-standard setup required by a derived class.
   */
  virtual void perform_nonstandard_setup() {}
#ifdef IS_PARALLEL
  /** \brief MPI communicator */
  MPI_Comm mpi_communicator;

  /** \brief Locally owned DoF index set */
  IndexSet locally_owned_dofs;

  /** \brief Locally relevant DoF index set */
  IndexSet locally_relevant_dofs;
#endif

  /** \brief Conditional output stream, necessary for parallel computations */
  ConditionalOStream cout1;
  ConditionalOStream cout2;

  /** \brief Timer */
  TimerOutput timer;

  /** \brief triangulation */
  LocalTriangulation triangulation;

  /** \brief name of problem */
  std::string problem_name;

  /** \brief input parameters for conservation law */
  const ConservationLawParameters<dim> parameters;

  /** \brief number of components in the system */
  const unsigned int n_components;

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
  LocalVector new_solution;
  /** \brief solution of previous time step */
  LocalVector old_solution;

  /** \brief solution step in Newton loop; vector for temporary storage */
  LocalVector solution_step;
  /** \brief new time (end of time step) */
  double new_time;
  /** \brief old time (beginning of time step) */
  double old_time;
  /** \brief system right-hand side */
  LocalVector system_rhs;

  /** \brief constrained sparsity pattern */
  SparsityPattern constrained_sparsity_pattern;
  /** \brief unconstrained sparsity pattern */
  SparsityPattern unconstrained_sparsity_pattern;
  /** \brief consistent mass matrix */
  LocalMatrix consistent_mass_matrix;
  /** \brief lumped mass matrix */
  LocalMatrix lumped_mass_matrix;
  /** \brief mass matrix to be used */
  const LocalMatrix * mass_matrix;
  /** \brief system matrix */
  LocalMatrix system_matrix;

  /** \brief number of Dirichlet boundaries */
  unsigned int n_dirichlet_boundaries;
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
  CellMap cell_diameter;
  CellMap max_flux_speed_cell;

  // viscosity
  std::shared_ptr<Viscosity<dim>> viscosity;
  std::shared_ptr<Viscosity<dim>> constant_viscosity;
  std::shared_ptr<Viscosity<dim>> low_order_viscosity;
  std::shared_ptr<Viscosity<dim>> entropy_viscosity;

  std::unique_ptr<ArtificialDiffusion<dim>> artificial_diffusion;

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
    std::vector<LocalVector> a;
    std::vector<double> b;
    std::vector<double> c;

    bool solution_computed_in_last_stage;
    bool is_explicit;

    std::vector<LocalVector> f;
  };

  /** \brief Runge-Kutta parameters */
  RungeKuttaParameters rk;

  /** \brief Default end time for test problem */
  double default_end_time;
  /** \brief Flag to signal that test problem has a default end time */
  bool has_default_end_time;
  /** \brief Chosen end time */
  double end_time;
};

#include "src/base/ConservationLaw.cc"

#endif
