#ifndef ScalarConservationLawFlux_h
#define ScalarConservationLawFlux_h

#include <deal.II/base/tensor.h>

using namespace dealii;

/**
 * \brief Abstract base class for scalar conservation law fluxes.
 */
template <int dim>
class ScalarConservationLawFlux
{
public:
  /**
   * \brief Constructor.
   *
   * \param[in] is_linear_  flag that flux is linear
   */
  ScalarConservationLawFlux(const bool & is_linear_) : is_linear(is_linear_) {}
  /** \brief Evaluates conservation law flux at a number of points. */
  virtual std::vector<Tensor<1, dim>> evaluate(
    const std::vector<double> & U, const std::vector<double> & x) const = 0;

  /** \brief Evaluates conservation law flux derivative at a number of points */
  virtual std::vector<Tensor<1, dim>> evaluate_derivative(
    const std::vector<double> & U, const std::vector<double> & x) const = 0;

  /** \brief flag that the flux is linear with respect to the solution. */
  const bool is_linear;
};

#endif
