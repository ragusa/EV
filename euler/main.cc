/**
 * \file main.cc
 * \brief The main source file.
 */

#include <iostream>
#include <deal.II/base/logstream.h>

#include "Burgers.h"
#include "BurgersParameters.h"
#include "Euler.h"
#include "EulerParameters.h"
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

    // declare parameter handler
    ParameterHandler parameter_handler;

#if defined(BURGERS)
    // read input and declare problem
    BurgersParameters<dimension>::declare_burgers_parameters(parameter_handler);
    parameter_handler.read_input("input_burgers");
    BurgersParameters<dimension> burgers_parameters;
    burgers_parameters.get_burgers_parameters(parameter_handler);

    // run problem
    Burgers<dimension> burgers_problem(burgers_parameters);
#elif defined(EULER)
    // read input and declare problem
    EulerParameters<dimension>::declare_parameters(parameter_handler);
    parameter_handler.read_input("input_euler");
    EulerParameters<dimension> parameters;
    parameters.get_parameters(parameter_handler);
    Euler<dimension> problem(parameters);
#elif defined(SHALLOWWATER)
    // ensure dimension <= 2
    Assert(dimension <= 2, ExcImpossibleInDim(dimension));

    // read input and declare problem
    ShallowWaterParameters<dimension>::declare_parameters(parameter_handler);
    parameter_handler.read_input("input_shallowwater");
    ShallowWaterParameters<dimension> parameters;
    parameters.get_parameters(parameter_handler);
    ShallowWater<dimension> problem(parameters);
#else
#error No valid conservation law defined in preprocessor
#endif

    // run problem
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
