#include <iostream>
#include "../../ExactSolutions.h"

using namespace dealii;

const unsigned int dim = 1; // number of spatial dimensions

int main(int argc, char ** argv) {
   try {

      // compute exact solution from function
      ThreeRegionExactSolution<dim> exact_solution;
      exact_solution.set_time(1.0);
      Point<dim> x(0.5);
      double function_value = exact_solution.value(x);

      // compute the expected value
      double incoming = 1.0;
      double sigma0 = 1.0;
      double sigma1 = 10.0;
      double q0 = 1.0;
      double q1 = 5.0;
      double s0 = 0.3;
      double s1 = 0.2;
      double ub = incoming*exp(-(sigma0*s0 + sigma1*s1));
      double uq0 = q0/sigma0*(1-exp(-sigma0*s0));
      double uq1 = q1/sigma1*(1-exp(-sigma1*s1));
      double uq = uq0*exp(-sigma1*s1) + uq1;
      double expected_value = ub + uq;

      // report the two values
      std::cout << "Computed value = " << function_value <<
         ", Expected value = " <<  expected_value << std::endl;
      std::cout << "Expected ub = " << ub << ", uq = " << uq << std::endl;

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
