#ifndef Executioner_cc
#define Executioner_cc

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/sparse_matrix.h>
#include "TransportParameters.h"

using namespace dealii;

/**
 * Class for executioner.
 */
template<int dim>
class Executioner
{
public:
  Executioner(const TransportParameters<dim> & parameters,
    const Triangulation<dim> & triangulation,
    const Tensor<1, dim> & transport_direction,
    const FunctionParser<dim> & cross_section_function,
    FunctionParser<dim> & source_function,
    Function<dim> & incoming_function);
  virtual ~Executioner();

  virtual void run() = 0;

  Vector<double> getFinalSolution() const;

protected:
  void setupSystem();
  void assembleInviscidSteadyStateMatrix();
  void assembleSteadyStateRHS(const double &t);

  const TransportParameters<dim> parameters;

  const FESystem<dim> fe;
  const FEValuesExtractors::Scalar flux;
  DoFHandler<dim> dof_handler;
  const unsigned int dofs_per_cell;
  const unsigned int n_cells;
  const QGauss<dim> cell_quadrature;
  const unsigned int n_q_points_cell;

  LinearSolver<dim> linear_solver;

  const Tensor<1, dim> transport_direction;
  const FunctionParser<dim> * const cross_section_function;
  FunctionParser<dim> * const source_function;
  Function<dim> * const incoming_function;

  ConstraintMatrix constraints;
  SparsityPattern constrained_sparsity_pattern;
  SparseMatrix<double> inviscid_ss_matrix;
  SparseMatrix<double> low_order_ss_matrix;
  SparseMatrix<double> low_order_diffusion_matrix;
  SparseMatrix<double> system_matrix;

  Vector<double> ss_rhs;
  Vector<double> new_solution;
};

#include "Executioner.cc"
#endif
