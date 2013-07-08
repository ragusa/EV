#include <iostream>

#include <deal.II/base/logstream.h>

#include "EulerEquationsBase.h"
#include "EulerEquationsBaseParameters.h"

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

      // get input parameters
      ParameterHandler                parameter_handler;
      EulerEquationsBaseParameters<1> parameters;

      parameters.declare_parameters(parameter_handler);
      parameter_handler.read_input(input_filename);
      parameters.get_parameters(parameter_handler);

      // run problem
      //EulerEquationsBase<1> euler_problem(parameters);
      //euler_problem.run();

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
