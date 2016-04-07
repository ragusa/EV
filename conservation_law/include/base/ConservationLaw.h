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
#include "include/diffusion/DomainInvariantDiffusion.h"
#include "include/diffusion/NoDiffusion.h"
#include "include/diffusion/GraphTheoreticDiffusion.h"
#include "include/diffusion/LaplacianDiffusion.h"
#include "include/entropy/Entropy.h"
#include "include/fct/ExplicitEulerFCT.h"
#include "include/fe/GradientMatrix.h"
#include "include/other/Exceptions.h"
#include "include/other/Utilities.h"
#include "include/parameters/RunParameters.h"
#include "include/parameters/ProblemParameters.h"
#include "include/postprocessing/PostProcessor.h"
#include "include/solvers/LinearSolver.h"
#include "include/time_integrators/SSPRKTimeIntegrator.h"
#include "include/viscosity/ConstantViscosity.h"
#include "include/viscosity/DomainInvariantViscosity.h"
#include "include/viscosity/DMPLowOrderViscosity.h"
#include "include/viscosity/EntropyViscosity.h"
#include "include/viscosity/HighOrderViscosity.h"
#include "include/viscosity/LowOrderViscosity.h"
#include "include/viscosity/MaxWaveSpeed.h"
#include "include/viscosity/Viscosity.h"
#include "include/viscosity/ViscosityMultiplier.h"
#include "include/riemann/StarState.h"

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
  ConservationLaw(const RunParameters & params,
                  const unsigned int & n_components,
                  const bool & is_linear);
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

  using TimeStepSizeOption = typename RunParameters::TimeStepSizeOption;
  using TemporalDiscretizationClassification =
    typename RunParameters::TemporalDiscretizationClassification;
  using SSPRKDiscretization = typename RunParameters::SSPRKDiscretization;
  using Scheme = typename RunParameters::Scheme;
  using LowOrderScheme = typename RunParameters::LowOrderScheme;
  using HighOrderScheme = typename RunParameters::HighOrderScheme;
  using ViscosityType = typename RunParameters::ViscosityType;
  using DiffusionType = typename RunParameters::DiffusionType;
  using LinearSolverType = typename RunParameters::LinearSolverType;

  void initialize_system();

  void setup_system();

  std::shared_ptr<ArtificialDiffusion<dim>> create_artificial_diffusion(
    const DiffusionType & diffusion_type);

  void compute_error_for_refinement();

  void refine_mesh();

  void update_cell_sizes();

  void assemble_mass_matrix();

  void solve_runge_kutta(PostProcessor<dim> & postprocessor);

  double compute_dt_from_cfl_condition();

  double compute_dt_from_dmp_cfl_condition();

  double compute_dt_from_di_cfl_condition();

  double compute_cfl_number(const double & dt) const;

  void check_nan();

  void output_viscosity(PostProcessor<dim> & postprocessor,
                        const bool & is_transient = false,
                        const double & time = 0.0);

  void perform_fct_ssprk_step(const double & dt,
                              const double & old_stage_dt,
                              const unsigned int & n,
                              const std::shared_ptr<ExplicitEulerFCT<dim>> & fct,
                              SSPRKTimeIntegrator<dim> & ssprk);

  void get_dirichlet_dof_indices(
    std::vector<unsigned int> & dirichlet_dof_indices);

  virtual void compute_inviscid_ss_matrix(const Vector<double> & solution,
                                          SparseMatrix<double> & matrix);

  virtual void compute_ss_reaction(Vector<double> & ss_reaction);

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

  /**
   * \brief Computes the steady-state flux vector.
   *
   * This function computes the steady-state flux vector \f$\mathrm{f}^n\f$:
   * \f[
   *   \mathrm{f}_i^n = \int\limits_{S_i}
   *     \varphi_i(\mathbf{x})\nabla\cdot\mathbf{f}(\tilde{\mathbf{u}}^n)dV \,.
   * \f]
   *
   * \param[in] dt time step size \f$\Delta t\f$
   * \param[in] solution solution vector \f$\mathrm{U}^n\f$
   * \param[out] ss_flux steady-state flux vector \f$\mathrm{f}^n\f$
   */
  virtual void compute_ss_flux(const double & dt,
                               const Vector<double> & solution,
                               Vector<double> & ss_flux) = 0;

  /**
   * \brief Computes the steady-state right hand side vector.
   *
   * This function computes the steady-state flux vector \f$\mathrm{q}^n\f$:
   * \f[
   *   \mathrm{q}_i^n = \int\limits_{S_i}
   *     \varphi_i(\mathbf{x}) q(t^n,\mathbf{x})dV \,.
   * \f]
   *
   * \param[in] t time at which to evaluate residual \f$t^n\f$
   * \param[out] ss_rhs steady-state rhs vector \f$\mathrm{q}^n\f$
   */
  virtual void compute_ss_rhs(const double & t, Vector<double> & ss_rhs) = 0;

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

  virtual std::shared_ptr<StarState<dim>> create_star_state() const;

  virtual std::shared_ptr<ExplicitEulerFCT<dim>> create_fct() const = 0;

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

  /** \brief input parameters for conservation law */
  const RunParameters parameters;

  /** \brief Pointer to problem parameters */
  ProblemParameters<dim> * problem_base_parameters;

  /** \brief number of components in the system */
  const unsigned int n_components;

  /** \brief Flag to signal if star states arise in Riemann problem of
   *         conservaiton law */
  const bool are_star_states;

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

  /** \brief Pointer to gradient matrix */
  std::shared_ptr<GradientMatrix<dim>> gradient_matrix;

  /** \brief Pointer to star state */
  std::shared_ptr<StarState<dim>> star_state;

  /** \brief Pointer to linear solver */
  std::shared_ptr<LinearSolver<dim>> linear_solver;

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
  /** \brief solution of previous previous time step */
  LocalVector older_solution;
  /** \brief old stage solution */
  LocalVector old_stage_solution;

  /** \brief solution step in Newton loop; vector for temporary storage */
  LocalVector solution_step;
  /** \brief new time */
  double new_time;
  /** \brief old time */
  double old_time;

  /** \brief system right-hand side */
  LocalVector system_rhs;
  /** \brief steady-state flux vector */
  LocalVector ss_flux;
  /** \brief steady-state reaction vector */
  LocalVector ss_reaction;
  /** \brief steady-state right-hand side vector */
  LocalVector ss_rhs;
  /** \brief temporary vector */
  LocalVector tmp_vector;

  /** \brief constrained sparsity pattern */
  SparsityPattern constrained_sparsity_pattern;
  /** \brief unconstrained sparsity pattern */
  SparsityPattern unconstrained_sparsity_pattern;

  /** \brief system matrix */
  LocalMatrix system_matrix;
  /** \brief consistent mass matrix */
  LocalMatrix consistent_mass_matrix;
  /** \brief lumped mass matrix */
  LocalMatrix lumped_mass_matrix;
  /** \brief Low-order diffusion matrix \f$\mathrm{D}^{L,n}\f$ */
  LocalMatrix low_order_diffusion_matrix;
  /** \brief High-order diffusion matrix \f$\mathrm{D}^{H,n}\f$ */
  LocalMatrix high_order_diffusion_matrix;
  /** \brief Inviscid steady-state matrix \f$\mathbf{A}\f$ */
  LocalMatrix inviscid_ss_matrix;
  /** \brief Low-order steady-state matrix \f$\mathbf{A}^L\f$ */
  LocalMatrix low_order_ss_matrix;

  /** \brief Vector of component names */
  std::vector<std::string> component_names;
  /** \brief Vector of data component interpretations (scalar or vector) */
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretations;

  /** \brief Minimum cell diameter */
  double minimum_cell_diameter;
  /** \brief Maximum flux speed in domain */
  double max_flux_speed;

  // maps
  CellMap max_flux_speed_cell;

  // viscosity
  /** \brief low-order viscosity type */
  ViscosityType low_order_viscosity_type;
  /** \brief entropy viscosity type */
  ViscosityType entropy_viscosity_type;
  /** \brief high-order viscosity type */
  ViscosityType high_order_viscosity_type;
  /** \brief low-order viscosity \f$\nu^L\f$ */
  std::shared_ptr<Viscosity<dim>> low_order_viscosity;
  /** \brief entropy viscosity \f$\nu^\eta\f$ */
  std::shared_ptr<Viscosity<dim>> entropy_viscosity;
  /** \brief high-order viscosity \f$\nu^H\f$ */
  std::shared_ptr<Viscosity<dim>> high_order_viscosity;

  // diffusion
  /** \brief low-order diffusion type */
  DiffusionType low_order_diffusion_type;
  /** \brief high-order diffusion type */
  DiffusionType high_order_diffusion_type;
  /** \brief low-order diffusion */
  std::shared_ptr<ArtificialDiffusion<dim>> low_order_diffusion;
  /** \brief high-order diffusion */
  std::shared_ptr<ArtificialDiffusion<dim>> high_order_diffusion;

  /** \brief Flag that conservation law is scalar */
  const bool is_scalar;
  /** \brief Flag that conservation law is linear */
  const bool is_linear;

  /** \brief Flag that the inviscid steady-state matrix needs to be computed */
  bool need_to_compute_inviscid_ss_matrix;

  /** \brief Flag to signal being in last adaptive refinement cycle */
  bool in_final_cycle;
  /** \brief Estimation of error per cell for adaptive mesh refinement */
  Vector<float> estimated_error_per_cell;

  /** \brief Indices of DoFs subject to Dirichlet boundary conditions */
  std::vector<unsigned int> dirichlet_dof_indices;
};

#include "src/base/ConservationLaw.cc"

#endif
