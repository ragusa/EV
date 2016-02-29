#ifndef PostProcessor_cc
#define PostProcessor_cc

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <sys/stat.h>
#include "TransportParameters.h"

using namespace dealii;

/**
 * Class for outputting solutions and evaluating error and convergence.
 */
template <int dim>
class PostProcessor
{
public:
  /** \brief Alias for temporal discretization */
  using TemporalDiscretization =
    typename TransportParameters<dim>::TemporalDiscretization;

  /** \brief Alias for theta method */
  using ThetaMethod = typename TransportParameters<dim>::ThetaMethod;

  /** \brief Alias for SSPRK method */
  using SSPRKMethod = typename TransportParameters<dim>::SSPRKMethod;

  PostProcessor(const TransportParameters<dim> & parameters,
                const bool has_exact_solution,
                std::shared_ptr<Function<dim>> & exact_solution_function);
  ~PostProcessor();

  void output_results(const Vector<double> & solution,
                      const DoFHandler<dim> & dof_handler,
                      const Triangulation<dim> & triangulation);
  void output_solution(const Vector<double> & solution,
                       const DoFHandler<dim> & dof_handler,
                       const std::string & output_string) const;
  void output_viscosity(const Vector<double> & low_order_viscosity,
                        const Vector<double> & entropy_viscosity,
                        const Vector<double> & high_order_viscosity,
                        const DoFHandler<dim> & dof_handler) const;
  void evaluate_error(const Vector<double> & solution,
                      const DoFHandler<dim> & dof_handler,
                      const Triangulation<dim> & triangulation);
  void update_dt(const double & dt);
  void setCycle(const unsigned int & cycle);
  bool askIfLastCycle() const;

private:
  void output_grid(const Triangulation<dim> & triangulation) const;
  void create_directory(const std::string & dir) const;

  const TransportParameters<dim> parameters;

  ConvergenceTable convergence_table;

  const bool has_exact_solution;
  std::shared_ptr<Function<dim>> exact_solution_function;

  double dt_nominal;
  const bool is_steady_state;

  const FESystem<dim> fe;

  std::string output_dir;
  std::string appendage_string;
  std::string filename_exact;

  const QGauss<dim> cell_quadrature;

  unsigned int current_cycle;
  bool is_last_cycle;
};

#include "PostProcessor.cc"
#endif
