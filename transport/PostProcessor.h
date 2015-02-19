#ifndef PostProcessor_cc
#define PostProcessor_cc

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function_parser.h>
//#include <deal.II/base/table_handler.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

using namespace dealii;

/** \brief Class for outputting solutions and evaluating error and convergence.
 */
template<int dim>
class PostProcessor {
   public:
      PostProcessor(
         const bool               &output_mesh,
         const bool               &output_exact_solution,
         const bool               &save_convergence_results,
         const bool               &has_exact_solution,
         FunctionParser<dim> &exact_solution_function,
         const double             &time,
         const double             &dt_nominal,
         const bool               &is_steady_state,
         const unsigned int       &refinement_option,
         const unsigned int       &final_refinement_level,
         const FESystem<dim>      &fe,
         const unsigned int       &degree,
         const unsigned int       &scheme_option,
         const unsigned int       &problem_ID,
         const QGauss<dim>        &cell_quadrature);
      ~PostProcessor();

      void output_results(const Vector<double>     &solution,
                          const DoFHandler<dim>    &dof_handler,
                          const Triangulation<dim> &triangulation);
      void output_solution(const Vector<double>  &solution,
                           const DoFHandler<dim> &dof_handler,
                           const std::string     &output_string,
                           const bool            &append_viscosity) const;
      void output_viscosity(const Vector<double>  &low_order_viscosity,
                            const Vector<double>  &entropy_viscosity,
                            const Vector<double>  &high_order_viscosity,
                            const DoFHandler<dim> &dof_handler) const;
      void evaluate_error(const Vector<double>     &solution,
                          const DoFHandler<dim>    &dof_handler,
                          const Triangulation<dim> &triangulation,
                          const unsigned int       &cycle);
      void update_dt(const double &dt);

   private:
      void output_grid(const Triangulation<dim> &triangulation) const;

      ConvergenceTable convergence_table;

      const bool               output_mesh;
      const bool               output_exact_solution;
      const bool               save_convergence_results;

      const bool               has_exact_solution;
      FunctionParser<dim>      *exact_solution_function;

      const double             time;
      double                   dt_nominal;
      const bool               is_steady_state;

      const unsigned int       refinement_option;
      const unsigned int       final_refinement_level;

      const FESystem<dim>      *fe;
      const unsigned int       degree;

      const unsigned int       scheme_option;
      const unsigned int       problem_ID;

      const QGauss<dim>        cell_quadrature;

      std::string              viscosity_string;
};

#include "PostProcessor.cc"
#endif
