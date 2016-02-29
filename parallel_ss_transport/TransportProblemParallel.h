#ifndef TransportProblem_cc
#define TransportProblem_cc

#include <deal.II/lac/generic_linear_algebra.h>
namespace LA
{
using namespace dealii::LinearAlgebraPETSc;
}

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>           // query MPI process info
#include <deal.II/base/conditional_ostream.h> // parallel cout
#include <deal.II/base/index_set.h>           // subset of indices for process

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_tools.h> // distribute_sparsity_pattern
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>

#include <deal.II/distributed/tria.h>            // distributed triangulation
#include <deal.II/distributed/grid_refinement.h> // distributed grid refinement

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include "TransportParameters.h"

using namespace dealii;

template <int dim>
class TransportProblem
{
public:
  TransportProblem(const TransportParameters<dim> & parameters);
  ~TransportProblem();
  void run();

private:
  // functions
  void initialize_system();
  void setup_system();
  void assemble_system();
  void solve();
  void refine_grid();
  void set_boundary_indicators();
  void process_problem_ID();
  void output_solution();

  // MPI communicator
  MPI_Comm mpi_communicator;

  // input parameters
  const TransportParameters<dim> & parameters;

  // mesh and dof data
  parallel::distributed::Triangulation<dim> triangulation;
  DoFHandler<dim> dof_handler;
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;
  unsigned int n_cells;
  unsigned int n_dofs;
  const unsigned int degree;
  const FE_Q<dim> fe;
  const unsigned int dofs_per_cell;
  const unsigned int faces_per_cell;

  // quadrature data
  const QGauss<dim> cell_quadrature;
  const QGauss<dim - 1> face_quadrature;
  const unsigned int n_q_points_cell;
  const unsigned int n_q_points_face;

  // sparse matrices, sparsity patterns, and constraints
  ConstraintMatrix constraints;
  SparsityPattern constrained_sparsity_pattern;
  LA::MPI::SparseMatrix system_matrix;

  // vectors for solutions and right hand sides
  LA::MPI::Vector solution;
  LA::MPI::Vector system_rhs;

  // physics data
  std::map<std::string, double> function_parser_constants;
  Tensor<1, dim> transport_direction;
  FunctionParser<dim> exact_solution_function;
  FunctionParser<dim> source_function;
  FunctionParser<dim> cross_section_function;
  FunctionParser<dim> incoming_function;
  std::string exact_solution_string;
  std::string source_string;
  std::string cross_section_string;
  std::string incoming_string;
  bool has_exact_solution;
  double x_min;
  double x_max;
  double domain_volume;

  // parallel output
  ConditionalOStream pcout;

  // timer
  TimerOutput computing_timer;
};

#include "TransportProblemParallel.cc"
#endif
