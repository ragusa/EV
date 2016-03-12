/**
 * \file main.cc
 * \brief The main source file.
 */

#include <iostream>
#include <deal.II/base/logstream.h>

#include "include/other/CMakeVars.h"
//#include "include/base/Burgers.h"
//#include "include/base/Euler.h"
#include "include/base/ShallowWater.h"
#include "include/base/Transport.h"
#include "include/parameters/BurgersRunParameters.h"
#include "include/parameters/EulerRunParameters.h"
#include "include/parameters/ShallowWaterRunParameters.h"
#include "include/parameters/TransportRunParameters.h"

using namespace dealii;

int main(int argc, char * argv[])
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
    BurgersRunParameters::declare_parameters(parameter_handler);
    if (argc > 1)
      parameter_handler.read_input(argv[1]);
    else
      parameter_handler.read_input("input/burgers.prm");
    BurgersRunParameters parameters;
    parameters.get_parameters(parameter_handler);
    Burgers<dimension> problem(parameters);
#elif defined(EULER)
    // read input and declare problem
    EulerRunParameters::declare_parameters(parameter_handler);
    if (argc > 1)
      parameter_handler.read_input(argv[1]);
    else
      parameter_handler.read_input("input/euler.prm");
    EulerRunParameters parameters;
    parameters.get_parameters(parameter_handler);
    Euler<dimension> problem(parameters);
#elif defined(SHALLOWWATER)
    // ensure dimension <= 2
    Assert(dimension <= 2, ExcImpossibleInDim(dimension));

    // read input and declare problem
    ShallowWaterRunParameters::declare_parameters(parameter_handler);
    if (argc > 1)
      parameter_handler.read_input(argv[1]);
    else
      parameter_handler.read_input("input/shallowwater.prm");
    ShallowWaterRunParameters parameters;
    parameters.get_parameters(parameter_handler);
    ShallowWater<dimension> problem(parameters);
#elif defined(TRANSPORT)
    // read input and declare problem
    TransportRunParameters::declare_parameters(parameter_handler);
    if (argc > 1)
      parameter_handler.read_input(argv[1]);
    else
      parameter_handler.read_input("input/transport.prm");
    TransportRunParameters parameters;
    parameters.get_parameters(parameter_handler);
    Transport<dimension> problem(parameters);
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
