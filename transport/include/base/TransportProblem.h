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
#include "CMakeVars.h"
#include "Exceptions.h"
#include "TransportParameters.h"
#include "TransportProblemParameters.h"
#include "PostProcessor.h"
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
  /** \brief Alias for temporal discretization */
  using TemporalDiscretization =
    typename TransportParameters<dim>::TemporalDiscretization;

  TransportProblem(const TransportParameters<dim> & parameters);

  void run();

private:
  void initializeSystem();

  void get_problem_parameters();

  /** \brief Conditional output stream 1 */
  ConditionalOStream cout1;

  /** \brief Conditional output stream 2 */
  ConditionalOStream cout2;

  /** \brief run parameters */
  const TransportParameters<dim> parameters;

  /** \brief flag that problem is time-dependent */
  const bool is_time_dependent;

  /** \brief problem parameters */
  TransportProblemParameters<dim> problem_parameters;

  /** \brief timer */
  TimerOutput timer;

  /** \brief mesh */
  Triangulation<dim> triangulation;
};

#include "TransportProblem.cc"
#endif
