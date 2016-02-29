/**
 * \file BurgersMaxWaveSpeed.h
 * \brief Provides the header for the BurgersMaxWaveSpeed class.
 */

#ifndef BurgersMaxWaveSpeed_cc
#define BurgersMaxWaveSpeed_cc

#include "include/viscosity/MaxWaveSpeed.h"

using namespace dealii;

/**
 * \brief Class for computing the maximum wave speed for the Burgers
 *        equation, used in the computation of the domain-invariant viscosity.
 */
template <int dim>
class BurgersMaxWaveSpeed : public MaxWaveSpeed<dim>
{
public:
  BurgersMaxWaveSpeed();

  double compute(const std::vector<double> & solution_left,
                 const std::vector<double> & solution_right,
                 const Tensor<1, dim> & normal_vector) const override;
};

#include "src/viscosity/BurgersMaxWaveSpeed.cc"

#endif
