#ifndef TransportProblem_cc
#define TransportProblem_cc

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
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
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include "TransportParameters.h"
#include "PostProcessor.h"
#include "ExactSolutions.h"
#include "RefinementHandler.h"
#include "SteadyStateExecutioner.h"
#include "TransientExecutioner.h"

using namespace dealii;

/**
 * Class for defining and running a transport problem.
 */
template <int dim>
class TransportProblem
{
public:
  TransportProblem(const TransportParameters<dim> & parameters);
  ~TransportProblem();

  void run();

private:
  void initializeSystem();
  void processProblemID();

  // input parameters
  const TransportParameters<dim> parameters;

  // mesh
  Triangulation<dim> triangulation;

  // physics data
  const bool is_time_dependent;
  std::map<std::string, double> function_parser_constants;
  Tensor<1, dim> transport_direction;
  FunctionParser<dim> initial_conditions;
  std::shared_ptr<Function<dim>> exact_solution_function;
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

  // timer
  TimerOutput timer;
};

#include "TransportProblem.cc"
#endif
