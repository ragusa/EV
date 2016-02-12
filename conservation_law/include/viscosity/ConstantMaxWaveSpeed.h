/**
 * \file ConstantMaxWaveSpeed.h
 * \brief Provides the header for the ConstantMaxWaveSpeed class.
 */

#ifndef ConstantMaxWaveSpeed_cc
#define ConstantMaxWaveSpeed_cc

#include "include/viscosity/MaxWaveSpeed.h"

using namespace dealii;

/**
 * \brief Class for computing the maximum wave speed that is constant.
 */
template <int dim>
class ConstantMaxWaveSpeed : public MaxWaveSpeed<dim>
{
public:
  ConstantMaxWaveSpeed(const double & constant_speed);

  double compute(const std::vector<double> & solution_left,
                 const std::vector<double> & solution_right,
                 const Tensor<1, dim> & normal_vector) const override;

private:
  const double constant_speed;
};

#include "src/viscosity/ConstantMaxWaveSpeed.cc"

#endif
