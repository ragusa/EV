/**
 * \file ShallowWaterStarState.h
 * \brief Provides the header for the ShallowWaterStarState class.
 */

#ifndef ShallowWaterStarState_cc
#define ShallowWaterStarState_cc

#include "include/viscosity/StarState.h"

using namespace dealii;

/**
 * \brief Class for computing the star state for the shallow water equations
 *
 * The velocity in the star region is computed as
 * \f[
 *   u^* = \frac{1}{2}\left(u_L + u_R + \mathcal{W}_R(h^*,\mathbf{u}_R)
 *     - \mathcal{W}_R(h^*,\mathbf{u}_R)\right)
 * \f]
 */
template <int dim>
class ShallowWaterStarState : public StarState<dim>
{
public:
  ShallowWaterStarState(const double & gravity);

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

#include "src/viscosity/ShallowWaterStarState.cc"

#endif
