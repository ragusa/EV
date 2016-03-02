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

#include "EntropyViscosity.h"
#include "FCT.h"
#include "NonlinearSolver.h"
#include "PostProcessor.h"
#include "TransportRunParameters.h"

using namespace dealii;

/**
 * \brief Abstract base class for executioner of a transport problem.
 */
template <int dim>
class TransportExecutioner
{
public:
  TransportExecutioner(const TransportRunParameters<dim> & parameters,
              Triangulation<dim> & triangulation,
              const Tensor<1, dim> & transport_direction,
              const double & transport_speed,
              const FunctionParser<dim> & cross_section_function,
              FunctionParser<dim> & source_function,
              Function<dim> & incoming_function,
              const double & domain_volume,
              PostProcessor<dim> & postprocessor);

  virtual void run() = 0;

  Vector<double> getFinalSolution() const;

  void print_solution() const;

protected:
  void setupSystem();
  void assembleInviscidSteadyStateMatrix();
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

  const TransportRunParameters<dim> parameters;

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

  LinearSolver<dim> linear_solver;

  NonlinearSolver<dim> nonlinear_solver;

  const Tensor<1, dim> transport_direction;
  const double transport_speed;
  const FunctionParser<dim> * const cross_section_function;
  FunctionParser<dim> * const source_function;
  Function<dim> * const incoming_function;
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
  Vector<double> ss_rhs;
  Vector<double> new_solution;
  Vector<double> cumulative_antidiffusion;

  std::vector<unsigned int> dirichlet_nodes;

  PostProcessor<dim> * const postprocessor;
};

#include "TransportExecutioner.cc"
#endif
