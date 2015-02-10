#ifndef TransportProblem_cc
#define TransportProblem_cc

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/table.h>

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
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include "TransportParameters.h"
#include "EntropyViscosity.h"
#include "LinearSolver.h"
#include "SSPRungeKuttaTimeIntegrator.h"

using namespace dealii;

/** \brief Class for defining a transport problem
 */
template<int dim>
class TransportProblem {
   public:
      TransportProblem(const TransportParameters &parameters);
      ~TransportProblem();
      void run();

   private:
      // main functions
      void initialize_system();
      void setup_system();
      void assemble_mass_matrices();
      void assemble_inviscid_ss_matrix();
      void assemble_ss_rhs(const double &t);
      void solve_steady_state(const LinearSolver<dim> &linear_solver);
      void refine_grid();

      // post-processing
      void output_results();
      void output_solution(const Vector<double>  &solution,
                           const DoFHandler<dim> &dof_handler,
                           const std::string     &output_string,
                           const bool            &append_viscosity) const;
      void output_grid() const;
      void evaluate_error(const unsigned int cycle);

      // input parameters
      const TransportParameters &parameters;

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
      SparseMatrix<double> low_order_viscous_matrix;
      SparseMatrix<double> high_order_viscous_matrix;
      SparseMatrix<double> consistent_mass_matrix;
      SparseMatrix<double> lumped_mass_matrix;
//      SparseMatrix<double> flux_correction_matrix;

      // vectors for solutions and right hand sides
      Vector<double> new_solution;
      Vector<double> old_solution;
      Vector<double> older_solution;
      Vector<double> system_rhs;
      Vector<double> ss_rhs;

      // viscosity vectors
      Vector<double> entropy_viscosity;
      Vector<double> low_order_viscosity;
      Vector<double> high_order_viscosity;

      // convergence table
      ConvergenceTable convergence_table;

      // physics data
      void process_problem_ID();
      std::map<std::string,double> function_parser_constants;
      Tensor<1,dim> transport_direction;
      FunctionParser<dim> initial_conditions;
      FunctionParser<dim> exact_solution;
      FunctionParser<dim> source_function;
      FunctionParser<dim> cross_section_function;
      FunctionParser<dim> incoming_function;
      std::string initial_conditions_string;
      std::string exact_solution_string;
      std::string source_string;
      std::string cross_section_string;
      std::string incoming_string;
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

      // low-order max-principle viscosity functions and data
      void compute_viscous_bilinear_forms();
      void compute_low_order_viscosity();
      void compute_viscous_matrix(const Vector<double> &viscosity,
                                  SparseMatrix<double> &viscous_matrix);
      SparsityPattern unconstrained_sparsity_pattern;
      SparseMatrix<double> viscous_bilinear_forms;

/*
      // high-order max-principle viscosity functions and data
      void assemble_flux_correction_matrix(const double &dt);
      void compute_limiting_coefficients();
      void get_matrix_row(const SparseMatrix<double>      &matrix,
                          const unsigned int              &i,
                                std::vector<double>       &row_values,
                                std::vector<unsigned int> &row_indices,
                                unsigned int              &n_col
                         );
      Vector<double> R_plus;
      Vector<double> R_minus;

      // FCT quantities
      Vector<double> flux_correction_vector;
      Vector<double> min_values;
      Vector<double> max_values;
      Vector<double> Q_plus;
      Vector<double> Q_minus;
      // max-principle checks
      bool check_max_principle(const double &dt,
                               const bool   &using_high_order);
      void debug_max_principle_low_order (const unsigned int &i, const double &dt);
      void compute_max_principle_bounds(const double &dt);
      void compute_steady_state_max_principle_bounds();
*/

      // Dirichlet boundary condition
      void get_dirichlet_nodes();
      std::vector<unsigned int> dirichlet_nodes;
};

#include "TransportProblem.cc"
#endif
