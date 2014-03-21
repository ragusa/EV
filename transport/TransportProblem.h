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
      // main functions
      void initialize_system();
      void process_problem_ID();
      void setup_system();
      void set_boundary_indicators();
      void assemble_mass_matrices();
      void assemble_system();
      void apply_Dirichlet_BC();
      void solve_step();
      void solve_linear_system(const SparseMatrix<double> &A,
                               const Vector<double> &b);
      void refine_grid();
      void output_results();
      void output_grid() const;
      void evaluate_error(const unsigned int cycle);
      double compute_viscosity(const typename DoFHandler<dim>::active_cell_iterator &cell,
                               const unsigned int  &i_cell,
                               const FEValues<dim> &fe_values,
                               const std::vector<double> &total_cross_section);

      // low-order max-principle viscosity functions and data
      void compute_viscous_bilinear_forms();
      void add_max_principle_viscous_bilinear_form();
      void compute_max_principle_viscosity();
      SparsityPattern unconstrained_sparsity_pattern;
      SparseMatrix<double> max_principle_viscosity_numerators;
      SparseMatrix<double> viscous_bilinear_forms;

      // high-order max-principle viscosity functions and data
      void assemble_high_order_quantities();
      void compute_limiting_coefficients();
      void compute_high_order_solution();
      void get_matrix_row(const SparseMatrix<double>      &matrix,
                          const unsigned int              &i,
                                std::vector<double>       &row_values,
                                std::vector<unsigned int> &row_indices,
                                unsigned int              &n_col
                         );
      SparseMatrix<double> auxiliary_mass_matrix;         // B matrix
      SparseMatrix<double> high_order_coefficient_matrix; // A matrix
      SparseMatrix<double> limiting_coefficient_matrix;   // L matrix
      Vector<double> high_order_viscosity;

      // input parameters
      const TransportParameters &parameters;

      // mesh and dof data
      Triangulation<dim> triangulation;
      DoFHandler<dim> dof_handler;
      const unsigned int degree;
      const FESystem<dim> fe;
      const FEValuesExtractors::Scalar flux;
      const unsigned int dofs_per_cell;

      // quadrature data
      const QGauss<dim>   cell_quadrature;
      const QGauss<dim-1> face_quadrature;
      const unsigned int n_q_points_cell;
      const unsigned int n_q_points_face;

      // sparse matrices, sparsity patterns, and constraints
      ConstraintMatrix constraints;
      SparsityPattern constrained_sparsity_pattern;
      SparseMatrix<double> system_matrix;
      SparseMatrix<double> consistent_mass_matrix;
      SparseMatrix<double> lumped_mass_matrix;

      // vectors for solutions and right hand sides
      Vector<double> old_solution;
      Vector<double> new_solution;
      Vector<double> system_rhs;
      Vector<double> ss_rhs;

      Vector<double> entropy_viscosity;
      Vector<double> low_order_viscosity;

      unsigned int nonlinear_iteration;

      ConvergenceTable convergence_table;

      // physics data
      Tensor<1,dim> transport_direction;
      const unsigned int incoming_boundary;
      FunctionParser<dim> initial_conditions;
      std::string initial_conditions_string;
      bool has_exact_solution;
      FunctionParser<dim> exact_solution;
      std::string exact_solution_string;
      unsigned int source_option;
      double source_value;
      unsigned int cross_section_option;
      double cross_section_value;
      double incoming_flux_value;
      bool is_linear;

      // entropy viscosity functions and data
      void compute_entropy_domain_average();
      double domain_volume;
      double domain_averaged_entropy;
      double max_entropy_deviation_domain;

      // CFL condition functions and data
      void check_CFL_condition(const double &dt) const;
      double minimum_cell_diameter;
   
      // checks
      void check_solution_nonnegative() const;
      bool check_local_discrete_max_principle() const;
};

#include "TransportProblem.cc"
#endif
