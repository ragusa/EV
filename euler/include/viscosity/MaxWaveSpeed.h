/**
 * \file MaxWaveSpeed.h
 * \brief Provides the header for the MaxWaveSpeed class.
 */

#ifndef MaxWaveSpeed_cc
#define MaxWaveSpeed_cc

using namespace dealii;

/**
 * \brief Abstract base class for computing the maximum wave speed, used in
 *        the computation of domain-invariant viscosity.
 */
template <int dim>
class MaxWaveSpeed
{
public:
  MaxWaveSpeed();

  virtual double compute(const std::vector<double> & solution_left,
                         const std::vector<double> & solution_right,
                         const Tensor<1, dim> & normal_vector) const = 0;
};

#endif
