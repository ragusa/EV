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
#include <cstdio>
#include <dirent.h>
#include <regex>
#include <sys/stat.h>
#include "ConservationLawParameters.h"
#include "Exceptions.h"

using namespace dealii;

/**
 * \brief Class for outputting solutions and evaluating error and convergence.
 */
template <int dim>
class PostProcessor
{
public:
  /** \brief Typedef for cell iterators */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  /** \brief Typedef for cell iterator map to double */
  typedef std::map<Cell, double> CellMap;

  PostProcessor(
    const ConservationLawParameters<dim> & parameters,
    const double & end_time,
    const bool has_exact_solution,
    std::shared_ptr<Function<dim>> & exact_solution_function,
    const std::string & problem_name,
    const std::vector<std::string> & solution_component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      solution_component_interpretations,
    const Triangulation<dim> & triangulation);

  ~PostProcessor();

  void output_results(
    const Vector<double> & solution,
    const DoFHandler<dim> & dof_handler,
    const Triangulation<dim> & triangulation,
    const std::shared_ptr<DataPostprocessor<dim>> data_postprocessor = nullptr);

  void output_solution(
    const Vector<double> & solution,
    const double & time,
    const DoFHandler<dim> & dof_handler,
    const std::string & output_string,
    const bool & output_1d_vtu = false,
    const std::shared_ptr<DataPostprocessor<dim>> data_postprocessor = nullptr);

  void output_solution_transient(
    const Vector<double> & solution,
    const double & time,
    const DoFHandler<dim> & dof_handler,
    const std::string & output_string,
    const bool & force_output = false,
    const std::shared_ptr<DataPostprocessor<dim>> aux_postprocessor = nullptr);

  void evaluate_error(const Vector<double> & solution,
                      const DoFHandler<dim> & dof_handler,
                      const Triangulation<dim> & triangulation);

  void update_dt(const double & dt);

  void set_cycle(const unsigned int & cycle);

  void output_cell_map(const CellMap & cell_map,
                       const double & time,
                       const std::string & quantity_string,
                       const DoFHandler<dim> & dof_handler) const;

private:
  void output_exact_solution(const double & time,
                             const std::shared_ptr<DataPostprocessor<dim>>
                               aux_postprocessor = nullptr);

  void output_at_dof_points(
    const Vector<double> & values,
    const double & time,
    const std::vector<std::string> & component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      component_interpretations,
    const DoFHandler<dim> & dof_handler,
    const std::string & output_string,
    const bool & output_1d_vtu = false,
    const std::shared_ptr<DataPostprocessor<dim>> data_postprocessor = nullptr,
    const std::string & transient_appendage = "");

  void output_convergence_data();

  void output_function(
    Function<dim> & function,
    const double & time,
    const std::vector<std::string> & component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      component_interpretations,
    const std::string & filename,
    const std::shared_ptr<DataPostprocessor<dim>> aux_postprocessor = nullptr);

  void output_grid(const Triangulation<dim> & triangulation) const;

  void create_directory(const std::string & dir) const;

  void remove_vtu_files(const std::string & directory) const;

  const ConservationLawParameters<dim> parameters;

  const double end_time;

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
  std::string timedisc_string;
  std::string appendage_string;
  std::string filename_exact;

  const QGauss<dim> cell_quadrature;

  unsigned int current_cycle;
  bool is_last_cycle;

  Triangulation<dim> fine_triangulation;
  DoFHandler<dim> fine_dof_handler;
  void createFineTriangulationAndDoFHandler(
    const Triangulation<dim> & triangulation);

  /** \brief Estimate of size of all transient output files */
  unsigned int transient_output_size;

  /** \brief Number used in the next transient solution file name */
  unsigned int transient_file_number;

  /** \brief Counter for the transient solution, i.e., the time step index,
   *         used in determining if a transient solution will be output */
  unsigned int transient_counter;

  /**
   * \brief Flag to signal that the transient solution was not output this
   *        time step.
   *
   * This is used because the final solution should be output in the transient,
   * even if the time step number was not scheduled to be output by the
   * user-specified output frequency.
   */
  bool transient_solution_not_output_this_step;

  /** \brief vector of times and corresponding file names */
  std::vector<std::pair<double, std::string>> times_and_filenames;
};

#include "PostProcessor.cc"
#endif
