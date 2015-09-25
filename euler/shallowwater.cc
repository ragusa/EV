/**
 * \file shallowwater.cc
 * \brief Provides the driver for solving the shallow water equations.
 */

#include <iostream>
#include <deal.II/base/logstream.h>
#include "ShallowWater.h"
#include "ShallowWaterParameters.h"

using namespace dealii;

int main(int, char **)
{
  try
  {
    dealii::deallog.depth_console(0);

    // spatial dimensions
    const int dimension = 1;

    // ensure dimension <= 2
    Assert(dimension <= 2, ExcImpossibleInDim(dimension));

    // declare input parameters and read them from input file
    ParameterHandler parameter_handler;
    ShallowWaterParameters<dimension>::declare_parameters(parameter_handler);
    parameter_handler.read_input("input_shallowwater");
    ShallowWaterParameters<dimension> parameters;
    parameters.get_parameters(parameter_handler);

    // run problem
    ShallowWater<dimension> problem(parameters);
    problem.run();
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
