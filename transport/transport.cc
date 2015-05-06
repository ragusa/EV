/** \brief solves the radiative transfer equation:
           \f[
              \frac{\partial\psi}{\partial t} 
              + \mathbf{\omega}\cdot\nabla\psi
              + \sigma_t \psi
              = q
           \f]
 */

#include <sstream>
#include <cstdlib>
#include <deal.II/base/parameter_handler.h>
#include "TransportParameters.h"
#include "TransportProblem.h"

using namespace dealii;

const unsigned int dim = 1; // number of spatial dimensions

/** \brief reads input file and then runs problem
 */
int main(int argc, char ** argv) {
   try {
      deallog.depth_console(0);

      // get input parameters
      ParameterHandler parameter_handler;
      TransportParameters<dim> parameters;
      parameters.declare_parameters(parameter_handler);
      parameter_handler.read_input("input");
      parameters.get_parameters(parameter_handler);

      // run problem
      TransportProblem<dim> transport_problem(parameters);
      transport_problem.run();

   } catch (std::exception &exc) {
      std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what()
            << std::endl << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;
      return 1;
   } catch (...) {
      std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;
      return 1;
   };

   return 0;
}
