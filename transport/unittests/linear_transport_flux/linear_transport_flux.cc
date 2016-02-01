#include <iostream>
#include <vector>
#include <memory>
#include <deal.II/base/tensor.h>
#include "../../ConservationLawFlux.h"
#include "../../LinearTransportFlux.h"

using namespace dealii;

const unsigned int dim = 3; // number of spatial dimensions

int main(int argc, char ** argv)
{
  try
  {
    // create shared pointer for conservation law flux
    std::shared_ptr<ConservationLawFlux<dim>> flux;

    // create linear transport flux
    Tensor<1, dim> velocity;
    velocity[0] = 1.0;
    velocity[1] = 5.0;
    velocity[2] = 10.0;
    std::shared_ptr<LinearTransportFlux<dim>> linear_flux =
      std::make_shared<LinearTransportFlux<dim>>(velocity);

    // point base class shared pointer to derived class object
    flux = linear_flux;

    // create test input values
    unsigned int n = 5;
    std::vector<double> U = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> x = {1.0, 4.0, 9.0, 16.0, 25.0};

    // evaluate flux and derivative
    std::vector<Tensor<1, dim>> f_values = flux->evaluate(U, x);
    std::vector<Tensor<1, dim>> dfdu_values = flux->evaluate_derivative(U, x);

    // print computed flux values
    std::cout << "Computed flux values:" << std::endl;
    for (unsigned int i = 0; i < n; ++i)
    {
      std::cout << "i = " << i << ":";
      for (unsigned int d = 0; d < dim; ++d)
        std::cout << " " << f_values[i][d];
      std::cout << std::endl;
    }

    std::cout << std::endl;

    // print expected flux values
    std::cout << "Expected flux values:" << std::endl;
    for (unsigned int i = 0; i < n; ++i)
    {
      std::cout << "i = " << i << ":";
      for (unsigned int d = 0; d < dim; ++d)
        std::cout << " " << velocity[d] * U[i];
      std::cout << std::endl;
    }

    std::cout << std::endl;

    // print computed flux derivative values
    std::cout << "Computed flux derivative values:" << std::endl;
    for (unsigned int i = 0; i < n; ++i)
    {
      std::cout << "i = " << i << ":";
      for (unsigned int d = 0; d < dim; ++d)
        std::cout << " " << dfdu_values[i][d];
      std::cout << std::endl;
    }

    std::cout << std::endl;

    // print expected flux derivative values
    std::cout << "Expected flux derivative values:" << std::endl;
    for (unsigned int i = 0; i < n; ++i)
    {
      std::cout << "i = " << i << ":";
      for (unsigned int d = 0; d < dim; ++d)
        std::cout << " " << velocity[d];
      std::cout << std::endl;
    }
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
