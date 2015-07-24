#ifndef TransportProblem_cc
#define TransportProblem_cc

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>
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
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include "TransportParameters.h"
#include "LowOrderViscosity.h"
#include "EntropyViscosity.h"
#include "LinearSolver.h"
#include "SSPRKTimeIntegrator.h"
#include "PostProcessor.h"
#include "FCT.h"
#include "ExactSolutions.h"

using namespace dealii;

/**
 * Class for defining a transport problem.
 */
template<int dim>
class TransportProblem {
   public:
      TransportProblem(const TransportParameters<dim> &parameters);
      ~TransportProblem();
      void run();

   private:
      // main functions
      void initialize_system();
      void setup_system();
      void assemble_mass_matrices();
      void assemble_inviscid_ss_matrix();
      void assemble_ss_rhs(const double &t);
      void solve_steady_state(LinearSolver<dim> &linear_solver);
      void refine_grid();

      // input parameters
      const TransportParameters<dim> &parameters;

      // mesh and dof data
      Triangulation<dim> triangulation;
      unsigned int n_cells;
      DoFHandler<dim> dof_handler;
      unsigned int n_dofs;
      const unsigned int degree;
      const FESystem<dim> fe;
      const FEValuesExtractors::Scalar flux;
      const unsigned int dofs_per_cell;
      const unsigned int faces_per_cell;
      void set_boundary_indicators();

      // quadrature data
      const QGauss<dim>   cell_quadrature;
      const QGauss<dim-1> face_quadrature;
      const unsigned int n_q_points_cell;
      const unsigned int n_q_points_face;

      // sparse matrices, sparsity patterns, and constraints
      ConstraintMatrix constraints;
      SparsityPattern constrained_sparsity_pattern;
      SparseMatrix<double> system_matrix;
      SparseMatrix<double> low_order_ss_matrix;
      SparseMatrix<double> high_order_ss_matrix;
      SparseMatrix<double> inviscid_ss_matrix;
      SparseMatrix<double> low_order_diffusion_matrix;
      SparseMatrix<double> high_order_diffusion_matrix;
      SparseMatrix<double> consistent_mass_matrix;
      SparseMatrix<double> lumped_mass_matrix;

      // vectors for solutions and right hand sides
      Vector<double> new_solution;
      Vector<double> old_solution;
      Vector<double> older_solution;
      Vector<double> oldest_solution;
      Vector<double> old_stage_solution;
      Vector<double> system_rhs;
      Vector<double> ss_rhs;

      // physics data
      void process_problem_ID();
      std::map<std::string,double> function_parser_constants;
      Tensor<1,dim> transport_direction;
      FunctionParser<dim> initial_conditions;
      std::shared_ptr<Function<dim> > exact_solution_function;
      FunctionParser<dim> source_function;
      FunctionParser<dim> cross_section_function;
      FunctionParser<dim> incoming_function;
      std::string initial_conditions_string;
      std::string exact_solution_string;
      std::string source_string;
      std::string cross_section_string;
      std::string incoming_string;
      ExactSolutionOption exact_solution_option;
      bool has_exact_solution;
      bool source_time_dependent;
      double x_min;
      double x_max;
      double domain_volume;

      // time functions and data
      double dt_nominal;
      double CFL_nominal;
      double enforce_CFL_condition(double &dt);
      double minimum_cell_diameter;

      // Dirichlet boundary condition
      void get_dirichlet_nodes();
      std::vector<unsigned int> dirichlet_nodes;

      // timer
      TimerOutput timer;
};

#include "TransportProblem.cc"
#endif
