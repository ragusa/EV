#ifndef ConservationLawFlux_h
#define ConservationLawFlux_h

#include <deal.II/base/tensor.h>

using namespace dealii;

/**
 * Conservation law flux class.
 */
template<int dim>
class ConservationLawFlux
{

public:

/** Constructor. */
ConservationLawFlux(const bool & is_linear) :
    is_linear(is_linear)
{
}

/** Evaluates conservation law flux at a number of points. */
virtual std::vector<Tensor<1,dim> > evaluate(
  const std::vector<double> & U,
  const std::vector<double> & x) const = 0;

/** Evaluates conservation law flux derivative at a number of points */
virtual std::vector<Tensor<1,dim> > evaluate_derivative(
  const std::vector<double> & U,
  const std::vector<double> & x) const = 0;

/** flag for the conservation law flux being linear with respect to the solution. */
const bool is_linear;

};

#endif
