#ifndef TransportProblem_cc
#define TransportProblem_cc

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
//#include <deal.II/base/parsed_function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/table.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
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

#include "TotalSource.h"
#include "TotalCrossSection.h"
#include "TransportParameters.h"

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
      void initialize_system();
      void process_problem_ID();
      void setup_system();
      void set_boundary_indicators();
      void assemble_mass_matrix();
      void assemble_system();
      void solve_step();
      void solve_linear_system(const SparseMatrix<double> &A,
                               const Vector<double> &b);
      void refine_grid();
      void output_results();
      void output_grid() const;
      void evaluate_error(const unsigned int cycle);
      void compute_viscous_bilinear_forms();
      void compute_max_principle_viscosity();
      void check_solution_nonnegative() const;
      bool check_local_discrete_max_principle() const;

      const TransportParameters &parameters; // input parameters
      unsigned int degree;

      Triangulation<dim> triangulation;
      DoFHandler<dim> dof_handler;
      FESystem<dim> fe;
      const unsigned int dofs_per_cell;

      QGauss<dim>   cell_quadrature;
      QGauss<dim-1> face_quadrature;
      const unsigned int n_q_points_cell;
      const unsigned int n_q_points_face;

      ConstraintMatrix constraints;

      SparsityPattern constrained_sparsity_pattern;
      SparseMatrix<double> system_matrix;
      SparseMatrix<double> mass_matrix;
      SparsityPattern unconstrained_sparsity_pattern;
      SparseMatrix<double> max_principle_viscosity_numerators;
      SparseMatrix<double> viscous_bilinear_forms;

      Vector<double> old_solution;
      Vector<double> new_solution;
      Vector<double> system_rhs;
      Vector<double> ss_rhs;

      Vector<double> old_first_order_viscosity;
      Vector<double> entropy_viscosity;
      Vector<double> max_principle_viscosity;

      unsigned int nonlinear_iteration;

      ConvergenceTable convergence_table;

      Tensor<1,dim> transport_direction;
      const unsigned int incoming_boundary;

      FunctionParser<dim> initial_conditions;
      std::string initial_conditions_string;
      FunctionParser<dim> exact_solution;
      std::string exact_solution_string;

      bool is_linear;
      bool has_exact_solution;

      unsigned int source_option;
      double source_value;
      unsigned int cross_section_option;
      double cross_section_value;
      double incoming_flux_value;
};

#include "TransportProblem.cc"
#endif
