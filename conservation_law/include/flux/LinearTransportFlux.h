#ifndef LinearTransportFlux_h
#define LinearTransportFlux_h

#include <vector>
#include <deal.II/base/tensor.h>
#include "ScalarConservationLawFlux.h"

using namespace dealii;

/**
 * Linear transport flux class.
 */
template <int dim>
class LinearTransportFlux : public ScalarConservationLawFlux<dim>
{
public:
  /** Constructor. */
  LinearTransportFlux(const Tensor<1, dim> & velocity);

  /** Evaluates conservation law flux at a number of points. */
  std::vector<Tensor<1, dim>> evaluate(
    const std::vector<double> & U, const std::vector<double> & x) const override;

  /** Evaluates conservation law flux derivative at a number of points. */
  std::vector<Tensor<1, dim>> evaluate_derivative(
    const std::vector<double> & U, const std::vector<double> & x) const override;

private:
  /** constant velocity tensor **/
  const Tensor<1, dim> velocity;
};

#include "LinearTransportFlux.cc"
#endif
