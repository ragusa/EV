#include <iostream>
#include "../../ExactSolutions.h"

using namespace dealii;

const unsigned int dim = 1; // number of spatial dimensions

int main(int argc, char ** argv)
{
  try
  {
    // create constructor arguments
    const std::vector<double> interface_positions({0.3, 0.6});
    const std::vector<double> region_sources({1.0, 5.0, 20.0});
    const std::vector<double> region_sigmas({1.0, 40.0, 20.0});
    const std::vector<double> direction({1.0, 0.0, 0.0});
    const double incoming = 1.0;

    // compute exact solution from function
    MultiRegionExactSolution<dim> exact_solution(
      interface_positions, region_sources, region_sigmas, direction, incoming);
    exact_solution.set_time(1.0);
    Point<dim> x(0.5);
    double function_value = exact_solution.value(x);

    // compute the expected value
    double sigma0 = 1.0;
    double sigma1 = 40.0;
    double q0 = 1.0;
    double q1 = 5.0;
    double s0 = 0.3;
    double s1 = 0.2;
    double ub = incoming * exp(-(sigma0 * s0 + sigma1 * s1));
    double uq0 = q0 / sigma0 * (1 - exp(-sigma0 * s0));
    double uq1 = q1 / sigma1 * (1 - exp(-sigma1 * s1));
    double uq = uq0 * exp(-sigma1 * s1) + uq1;
    double expected_value = ub + uq;

    // report the two values
    std::cout << "Computed value = " << function_value
              << ", Expected value = " << expected_value << std::endl;
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
