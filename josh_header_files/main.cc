#include <iostream>

#include "EulerEquationsBase.h"
#include "EulerEquationsBaseParameters.h"

/**
 * \fn    main
 * \brief reads input file and then runs problem
 */
int main(int argc, char ** argv) {
   try {
      dealii::deallog.depth_console(0);
/*
      // input filename
      std::string input_filename;
      if (argc < 2) {
         std::cout << "Need to provide an argument to specify the input file.\n"
            << "The argument '#' specifies the input file 'test_problem_#.in'"
            << std::endl;
         std::exit(1);
      } else {
         std::stringstream input_filename_ss;
         input_filename_ss << "input/test_problem_" << argv[1] << ".in";
         input_filename = input_filename_ss.str();
      }

      // get input parameters
      dealii::ParameterHandler     parameter_handler;
      EulerEquationsBaseParameters parameters;

      parameters.declare_parameters(parameter_handler);
      parameter_handler.read_input(input_filename);
      parameters.get_parameters(parameter_handler);

      // run problem
      EulerEquationsBase<dim> euler_problem(parameters);
      euler_problem.run();
	  */

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
