/**
 * \file burgers.cc
 * \brief Provides the driver for solving the Burgers equation.
 */

#include <iostream>
#include <deal.II/base/logstream.h>
#include "Burgers.h"
#include "BurgersParameters.h"

#ifdef IS_PARALLEL
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#endif

int main(int, char **)
{
  try
  {
    using namespace dealii;

#ifdef IS_PARALLEL
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
#endif

    dealii::deallog.depth_console(0);

    // spatial dimensions
    const int dimension = 1;

    // declare input parameters and read them from input file
    ParameterHandler parameter_handler;
    BurgersParameters<dimension>::declare_burgers_parameters(parameter_handler);
    parameter_handler.read_input("input_burgers");
    BurgersParameters<dimension> burgers_parameters;
    burgers_parameters.get_burgers_parameters(parameter_handler);

    // run problem
    Burgers<dimension> burgers_problem(burgers_parameters);
    burgers_problem.run();
  }
  catch (std::exception & exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
