#include <sstream>
#include <cstdlib>

#include <deal.II/base/parameter_handler.h>   // parameter handler
#include <deal.II/base/utilities.h>           // MPI query
#include <deal.II/base/logstream.h>           // deallog
#include <deal.II/base/conditional_ostream.h> // parallel cout

#include "TransportParameters.h"      // transport problem parameters class
#include "TransportProblemSerial.h"   // serial transport problem class

using namespace dealii;

const unsigned int dim = 2; // number of spatial dimensions

/** \brief reads input file and then runs problem
 */
int main(int argc, char ** argv) {
   try {
      // initialize and finalize MPI
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      // print number of MPI processes
      unsigned int n_processes = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      ConditionalOStream pcout(std::cout,
         (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));
      pcout << "Running with " << n_processes << " MPI processes" << std::endl;

      // check that input file was supplied as command-line argument
      if (argc != 2) {
         pcout << "Usage: transport <input_file>" << std::endl;
         std::exit(EXIT_FAILURE);
      }
      // get input file name
      std::string input_file = argv[1];

      // set log depth
      deallog.depth_console(0);

      // get input parameters
      ParameterHandler parameter_handler;
      TransportParameters<dim> parameters;
      parameters.declare_parameters(parameter_handler);
      parameter_handler.read_input(input_file);
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
