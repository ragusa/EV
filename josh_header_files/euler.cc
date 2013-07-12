#include <iostream>

#include <deal.II/base/logstream.h>

#include "EulerEquations.h"
#include "EulerEquationsParameters.h"

using namespace dealii;

/**
 * \fn    euler
 * \brief reads input file and then runs problem
 */
int main(int argc, char ** argv) {
   try {
      dealii::deallog.depth_console(0);
      // input filename
      std::string input_filename;
      if (argc < 2) {
         std::cout << "Need to provide an argument to specify the input file.\n"
            << "The argument '#' specifies the input file 'test_problem_#.in'"
            << std::endl;
         std::exit(1);
      } else {
         std::stringstream input_filename_ss;
         input_filename_ss << "input/input_" << argv[1] << ".in";
         input_filename = input_filename_ss.str();
         std::cout << "Using the input file: " << input_filename << std::endl;
      }

      const int dimension = 2;

      // get input parameters
      ParameterHandler                     parameter_handler;
      //EulerEquationsParameters<dimension>  euler_parameters;
      //ConservationLawParameters<dimension> conservation_law_parameters(dimension+2);

      //euler_parameters.declare_parameters(parameter_handler);
      //conservation_law_parameters.declare_parameters(parameter_handler);
      /* we can declare Euler equations parameters using a static function;
         the conservation law parameters equivalent declare_parameters()
         function cannot be static because it depends on a non-static
         member variable: the number of components n_components, which is
         passed through the ConservationLaw contructor. Therefore, the
         conservation law parameters declare_parameters() function is
         called in the ConservationLaw constructor, as well as the
         input file reading with the parameter handler
      */
      EulerEquationsParameters<dimension>::declare_parameters(parameter_handler);
      ConservationLawParameters<dimension>::declare_parameters(parameter_handler);
      parameter_handler.read_input(input_filename);

      //euler_parameters.get_parameters(parameter_handler);
      //conservation_law_parameters.get_parameters(parameter_handler);

      // run problem
      EulerEquations<dimension> euler_problem(parameter_handler);
      euler_problem.test_run();

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
