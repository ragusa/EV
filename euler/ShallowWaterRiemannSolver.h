/**
 * \file ShallowWaterRiemannSolver.h
 * \brief Provides the header for the ShallowWaterRiemannSolver class.
 */

#ifndef ShallowWaterRiemannSolver_cc
#define ShallowWaterRiemannSolver_cc

#include <deal.II/base/function.h>

using namespace dealii;

/**
 * \brief Class for computing the exact solution for a shallow water equations
 *        Riemann problem.
 */
template <int dim>
class ShallowWaterRiemannSolver : public Function<dim>
{
public:
  ShallowWaterRiemannSolver(const double & h_left,
                            const double & u_left,
                            const double & h_right,
                            const double & u_right,
                            const double & gravity,
                            const double & x_interface);

  double value(const Point<dim> & p, const unsigned int component) const override;

private:
  void evaluate_fluxes(const double & h,
                       const double & h_star,
                       const double & c_star,
                       double & f,
                       double & f_deriv) const;

  const double h_left;
  const double u_left;
  const double h_right;
  const double u_right;
  const double gravity;
  const double x_interface;
  const double c_left;
  const double c_right;
  double h_star;
  double u_star;
  double c_star;
};

#include "ShallowWaterRiemannSolver.cc"
#endif
