/**
 * \file ShallowWaterMaxWaveSpeed.h
 * \brief Provides the header for the ShallowWaterMaxWaveSpeed class.
 */

#ifndef ShallowWaterMaxWaveSpeed_cc
#define ShallowWaterMaxWaveSpeed_cc

#include "include/viscosity/MaxWaveSpeed.h"
#include "include/riemann/StarState.h"

using namespace dealii;

/**
 * \brief Class for computing the maximum wave speed for the shallow water
 *        equations, used in the computation of the domain-invariant viscosity.
 */
template <int dim>
class ShallowWaterMaxWaveSpeed : public MaxWaveSpeed<dim>
{
public:
  ShallowWaterMaxWaveSpeed(const std::shared_ptr<StarState<dim>> & star_state,
                           const double & gravity);

  double compute(const std::vector<double> & solution_left,
                 const std::vector<double> & solution_right,
                 const Tensor<1, dim> & normal_vector) const override;

protected:
  double compute_height_star(const double & h_left,
                             const double & u_left,
                             const double & h_right,
                             const double & u_right) const;

  void evaluate_fluxes(const double & h,
                       const double & h_K,
                       const double & c_K,
                       double & f,
                       double & f_deriv) const;

  const double gravity;
};

#include "src/viscosity/ShallowWaterMaxWaveSpeed.cc"

#endif
