/**
 * \file PostProcessor.h
 * \brief Provides the header for the Postprocessor class.
 */

#ifndef PostProcessor_cc
#define PostProcessor_cc

#include <cstdio>
#include <dirent.h>
#include <regex>
#include <sys/stat.h>
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
    const Triangulation<dim> & triangulation,
    const std::shared_ptr<DataPostprocessor<dim>> aux_postprocessor = nullptr);

  ~PostProcessor();

  void output_results(const Vector<double> & solution,
                      const DoFHandler<dim> & dof_handler,
                      const Triangulation<dim> & triangulation);

  void output_solution(const Vector<double> & solution,
                       const double & time,
                       const DoFHandler<dim> & dof_handler,
                       const std::string & output_string,
                       const bool & output_1d_vtu = false);

  void output_solution_transient(const Vector<double> & solution,
                                 const double & time,
                                 const DoFHandler<dim> & dof_handler,
                                 const std::string & output_string,
                                 const bool & force_output = false);

  void output_viscosity_transient(
    const std::vector<std::shared_ptr<Viscosity<dim>>> & viscosities,
    const std::vector<std::string> & names,
    const double & time,
    const DoFHandler<dim> & dof_handler,
    const bool & force_output = false);

  void evaluate_error(const Vector<double> & solution,
                      const DoFHandler<dim> & dof_handler,
                      const Triangulation<dim> & triangulation);

  void update_dt(const double & dt);

  void set_cycle(const unsigned int & cycle);

  void output_cell_maps(const std::vector<CellMap *> & cell_maps,
                        const std::vector<std::string> & names,
                        const std::string & filename_base,
                        const double & time,
                        const DoFHandler<dim> & dof_handler,
                        const bool & output_1d_vtu = false,
                        const std::string & transient_appendage = "");

  void output_viscosity(
    const std::vector<std::shared_ptr<Viscosity<dim>>> & viscosity,
    const std::vector<std::string> & names,
    const double & time,
    const DoFHandler<dim> & dof_handler,
    const bool & output_1d_vtu = false,
    const std::string & transient_appendage = "");

private:
  void output_exact_solution(const double & time);

  void output_at_dof_points(
    const Vector<double> & values,
    const double & time,
    const std::vector<std::string> & component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      component_interpretations,
    const DoFHandler<dim> & dof_handler,
    const std::string & output_string,
    const bool & output_1d_vtu = false,
    const bool & is_solution = true,
    const std::string & transient_appendage = "");

  void output_convergence_data();

  void output_function(
    Function<dim> & function,
    const double & time,
    const std::vector<std::string> & component_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
      component_interpretations,
    const std::string & filename,
    const bool & is_solution = true);

  void output_grid(const Triangulation<dim> & triangulation) const;

  void create_directory(const std::string & dir) const;

  void remove_vtu_files(const std::string & directory,
                        const std::string & filename_base) const;

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
  unsigned int transient_solution_file_number;

  /** \brief Number used in the next transient viscosity file name */
  unsigned int transient_viscosity_file_number;

  /** \brief Counter for the transient solution, i.e., the time step index,
   *         used in determining if a transient solution will be output */
  unsigned int transient_counter;

  /**
   * \brief Flag to signal that the transient solution and/or viscosity
   *        was not output this time step.
   *
   * This is used because the final solution and viscosity should be output
   * in the transient,
   * even if the time step number was not scheduled to be output by the
   * user-specified output frequency.
   */
  bool transient_not_output_this_step;

  /** \brief vector of times and corresponding solution file names */
  std::vector<std::pair<double, std::string>> times_and_solution_filenames;

  /** \brief vector of times and corresponding viscosity file names */
  std::vector<std::pair<double, std::string>> times_and_viscosity_filenames;

  /** \brief post-processor for derived quantities */
  std::shared_ptr<DataPostprocessor<dim>> solution_aux_postprocessor;
};

#include "PostProcessor.cc"
#endif
