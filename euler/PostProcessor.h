/**
 * \file PostProcessor.h
 * \brief Provides the header for the Postprocessor class.
 */

#ifndef PostProcessor_cc
#define PostProcessor_cc

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/vector_tools.h>
#include <sys/stat.h>
#include "ConservationLawParameters.h"

using namespace dealii;

/**
 * \brief Class for outputting solutions and evaluating error and convergence.
 */
template <int dim>
class PostProcessor
{
public:
  PostProcessor(
    const ConservationLawParameters<dim> & parameters,
    const bool has_exact_solution,
    std::shared_ptr<Function<dim>> & exact_solution_function,
    const std::string & problem_name,
    const std::vector<std::string> & solution_component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      solution_component_interpretations,
    const Triangulation<dim> & triangulation);

  ~PostProcessor();

  void output_results(const Vector<double> & solution,
                      const DoFHandler<dim> & dof_handler,
                      const Triangulation<dim> & triangulation);

  void output_results(const Vector<double> & solution,
                      const DoFHandler<dim> & dof_handler,
                      const Triangulation<dim> & triangulation,
                      const DataPostprocessor<dim> & data_postprocessor);

  void output_solution(const Vector<double> & solution,
                       const DoFHandler<dim> & dof_handler,
                       const std::string & output_string) const;

  void output_solution(const Vector<double> & solution,
                       const DoFHandler<dim> & dof_handler,
                       const std::string & output_string,
                       const DataPostprocessor<dim> & data_postprocessor) const;

  void evaluate_error(const Vector<double> & solution,
                      const DoFHandler<dim> & dof_handler,
                      const Triangulation<dim> & triangulation);

  void update_dt(const double & dt);

  void set_cycle(const unsigned int & cycle);

private:

  void output_viscosity(const Vector<double> & low_order_viscosity,
                        const Vector<double> & entropy_viscosity,
                        const Vector<double> & high_order_viscosity,
                        const DoFHandler<dim> & dof_handler) const;

  void output_exact_solution();

  void output_at_dof_points(
    const Vector<double> & values,
    const std::vector<std::string> & component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      component_interpretations,
    const DoFHandler<dim> & dof_handler,
    const std::string & output_string) const;

  void output_at_dof_points(
    const Vector<double> & values,
    const std::vector<std::string> & component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      component_interpretations,
    const DoFHandler<dim> & dof_handler,
    const std::string & output_string,
    const DataPostprocessor<dim> & data_postprocessor) const;

  void output_convergence_data();

  void output_function(
    const Function<dim> & function,
    const std::vector<std::string> & component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      component_interpretations,
    const std::string & filename) const;

  void output_function(
    Function<dim> & function,
    const std::vector<std::string> & component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      component_interpretations,
    const std::string & filename,
    const double & end_time) const;

  void output_grid(const Triangulation<dim> & triangulation) const;

  void create_directory(const std::string & dir) const;

  const ConservationLawParameters<dim> parameters;

  const std::string problem_name;

  ConvergenceTable convergence_table;

  const bool has_exact_solution;
  std::shared_ptr<Function<dim>> exact_solution_function;

  /** \brief List of names of each solution component */
  const std::vector<std::string> solution_component_names;
  /** \brief List of type (scalar or vector) of each solution component */
  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    solution_component_interpretations;

  double dt_nominal;
  const bool is_steady_state;

  const FESystem<dim> fe;

  std::string output_dir;
  std::string appendage_string;
  std::string filename_exact;

  const QGauss<dim> cell_quadrature;

  unsigned int current_cycle;
  bool is_last_cycle;

  Triangulation<dim> fine_triangulation;
  DoFHandler<dim> fine_dof_handler;
  void createFineTriangulationAndDoFHandler(
    const Triangulation<dim> & triangulation);
};

#include "PostProcessor.cc"
#endif
