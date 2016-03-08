/**
 * Solves a scalar transport equation of the form
 * \f[
 *   \frac{1}{c}\frac{\partial u}{\partial t}
 *   + \mathbf{\Omega}\cdot\nabla u
 *   + \sigma u
 *   = q
 * \f]
 */

#include <sstream>
#include <cstdlib>
#include <deal.II/base/parameter_handler.h>
#include "TransportRunParameters.h"
#include "TransportProblem.h"
#include "Exceptions.h"

using namespace dealii;

const unsigned int dim = 1; // number of spatial dimensions

/**
 * Reads input file and runs problem.
 */
int main(int argc, char * argv[])
{
  try
  {
    deallog.depth_console(0);

    // read input and declare problem
    ParameterHandler parameter_handler;
    TransportRunParameters<dim>::declare_parameters(parameter_handler);
    if (argc > 1)
      parameter_handler.read_input(argv[1]);
    else
      parameter_handler.read_input("../conservation_law/input/transport.prm");
    TransportRunParameters<dim> parameters;
    parameters.get_parameters(parameter_handler);

    // run problem
    TransportProblem<dim> transport_problem(parameters);
    transport_problem.run();
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
  };

  return 0;
}
