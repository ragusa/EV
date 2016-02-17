/**
 * \file MaxWaveSpeed.h
 * \brief Provides the header for the MaxWaveSpeed class.
 */

#ifndef MaxWaveSpeed_cc
#define MaxWaveSpeed_cc

#include "include/riemann/StarState.h"

using namespace dealii;

/**
 * \brief Abstract base class for computing the maximum wave speed, used in
 *        the computation of domain-invariant viscosity.
 */
template <int dim>
class MaxWaveSpeed
{
public:
  MaxWaveSpeed(const std::shared_ptr<StarState<dim>> & star_state = nullptr);

  virtual double compute(const std::vector<double> & solution_left,
                         const std::vector<double> & solution_right,
                         const Tensor<1, dim> & normal_vector) const = 0;

  void reset_max_wave_speed_domain();

  double get_max_wave_speed_domain() const;

protected:
  /** \brief Max wave speed found in domain */
  double max_wave_speed_domain;

  /** \brief Pointer to star state */
  std::shared_ptr<StarState<dim>> star_state;
};

#include "src/viscosity/MaxWaveSpeed.cc"

#endif
