/**
 * \file euler.cc
 * \brief Provides the driver for solving the Euler equations.
 */

#include <iostream>
#include <deal.II/base/logstream.h>
#include "Euler.h"
#include "EulerParameters.h"

using namespace dealii;

int main(int, char **)
{
  try
  {
    dealii::deallog.depth_console(0);

    // spatial dimensions
    const int dimension = 1;

    // declare input parameters and read them from input file
    ParameterHandler parameter_handler;
    EulerParameters<dimension>::declare_euler_parameters(parameter_handler);
    parameter_handler.read_input("input_euler");
    EulerParameters<dimension> euler_parameters;
    euler_parameters.get_euler_parameters(parameter_handler);

    // run problem
    Euler<dimension> euler_problem(euler_parameters);
    euler_problem.run();
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
