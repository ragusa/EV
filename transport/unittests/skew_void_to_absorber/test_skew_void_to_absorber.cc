#include <iostream>
#include "../../ExactSolutions.h"

using namespace dealii;

const unsigned int dim = 2; // number of spatial dimensions

int main(int argc, char ** argv) {
   try {

      // create necessary constructor arguments
      const std::vector<double> interface_positions = {0.5};
      const std::vector<double> region_sources = {0.0,0.0};
      const std::vector<double> region_sigmas = {0.0,10.0};
      const std::vector<double> direction({
         1.0/sqrt(2.0),
         1.0/sqrt(3.0),
         1.0/sqrt(6.0)
      });
      const double incoming = 1.0;

      // compute exact solution from function
      MultiRegionExactSolution<dim> exact_solution(
        interface_positions,
        region_sources,
        region_sigmas,
        direction,
        incoming
      );
      exact_solution.set_time(100.0);
      Point<dim> x(0.9,0.6);
      double function_value = exact_solution.value(x);

      // compute the expected value
      double sy = 0.1;
      double sx = sqrt(3.0/2.0)*sy;
      double sz = sqrt(3.0/6.0)*sy;
      double s = sqrt(pow(sx,2) + pow(sy,2) + pow(sz,2));
      double expected_value = exp(-10*s);

      // report the two values
      std::cout << "Computed value = " << function_value <<
         ", Expected value = " <<  expected_value << std::endl;

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
