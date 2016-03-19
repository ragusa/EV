#ifndef Executioner_cc
#define Executioner_cc

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/sparse_matrix.h>

#include "include/diffusion/ArtificialDiffusion.h"
#include "include/diffusion/NoDiffusion.h"
#include "include/diffusion/GraphTheoreticDiffusion.h"
#include "include/entropy/TransportEntropy.h"
#include "include/parameters/TransportProblemParameters.h"
#include "include/parameters/TransportRunParameters.h"
#include "include/postprocessing/PostProcessor.h"
#include "include/solvers/LinearSolver.h"
#include "include/solvers/NonlinearSolver.h"
#include "include/viscosity/ConstantViscosity.h"
#include "include/viscosity/DMPLowOrderViscosity.h"
#include "include/viscosity/EntropyViscosity.h"
#include "include/viscosity/HighOrderViscosity.h"

using namespace dealii;

/**
 * \brief Abstract base class for executioner of a transport problem.
 */
template <int dim>
class TransportExecutioner
{
public:
  using Scheme = typename RunParameters::Scheme;
  using LowOrderScheme = typename RunParameters::LowOrderScheme;
  using HighOrderScheme = typename RunParameters::HighOrderScheme;
  using ViscosityType = typename RunParameters::ViscosityType;
  using DiffusionType = typename RunParameters::DiffusionType;

  TransportExecutioner(const TransportRunParameters & parameters,
                       TransportProblemParameters<dim> & problem_parameters,
                       Triangulation<dim> & triangulation,
                       PostProcessor<dim> & postprocessor);

  virtual void run() = 0;

  Vector<double> getFinalSolution() const;

  void print_solution() const;

protected:
  void setupSystem();
  void assembleInviscidSteadyStateMatrix();
  void compute_reaction_vector();
  void assembleSteadyStateRHS(Vector<double> & rhs, const double & t);
  void setBoundaryIndicators();
  void getDirichletNodes();
  void applyDirichletBC(SparseMatrix<double> & A,
                        Vector<double> & b,
                        Vector<double> & x,
                        const double & t = 0.0);

  /** \brief Conditional output stream 1 */
  ConditionalOStream cout1;

  /** \brief Conditional output stream 2 */
  ConditionalOStream cout2;

  const TransportRunParameters parameters;

  const TransportProblemParameters<dim> * const problem_parameters;

  Triangulation<dim> * const triangulation;

  const FESystem<dim> fe;
  const FEValuesExtractors::Scalar flux;
  DoFHandler<dim> dof_handler;
  unsigned int n_dofs;
  const unsigned int dofs_per_cell;
  const unsigned int n_cells;
  const QGauss<dim> cell_quadrature;
  const QGauss<dim - 1> face_quadrature;
  const unsigned int n_q_points_cell;

  const Tensor<1, dim> transport_direction;
  const double transport_speed;
  const FunctionParser<dim> * const cross_section_function;
  FunctionParser<dim> * const source_function;
  std::shared_ptr<Function<dim>> incoming_function;
  const double domain_volume;

  ConstraintMatrix constraints;
  SparsityPattern constrained_sparsity_pattern;
  SparseMatrix<double> system_matrix;
  SparseMatrix<double> inviscid_ss_matrix;
  SparseMatrix<double> low_order_ss_matrix;
  SparseMatrix<double> high_order_ss_matrix;
  SparseMatrix<double> low_order_diffusion_matrix;
  SparseMatrix<double> high_order_diffusion_matrix;

  Vector<double> system_rhs;
  Vector<double> inviscid_ss_flux;
  Vector<double> reaction_vector;
  Vector<double> ss_rhs;
  Vector<double> new_solution;
  Vector<double> antidiffusion_vector;

  std::vector<unsigned int> dirichlet_nodes;

  LinearSolver<dim> linear_solver;

  NonlinearSolver<dim> nonlinear_solver;

  PostProcessor<dim> * const postprocessor;

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

  /** \brief low-order diffusion type */
  DiffusionType low_order_diffusion_type;
  /** \brief high-order diffusion type */
  DiffusionType high_order_diffusion_type;
  /** \brief low-order diffusion */
  std::shared_ptr<ArtificialDiffusion<dim>> low_order_diffusion;
  /** \brief high-order diffusion */
  std::shared_ptr<ArtificialDiffusion<dim>> high_order_diffusion;

private:
  std::shared_ptr<ArtificialDiffusion<dim>> create_artificial_diffusion(
    const DiffusionType & diffusion_type);
};

#include "src/executioners/TransportExecutioner.cc"
#endif
