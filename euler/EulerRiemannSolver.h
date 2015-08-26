#ifndef EulerRiemannSolver_cc
#define EulerRiemannSolver_cc

#include <deal.II/base/function.h>

using namespace dealii;

/**
 * Class for computing the exact solution for an Euler Riemann problem.
 */
template<int dim>
class EulerRiemannSolver : public Function<dim>
{
public:

  EulerRiemannSolver(
    const double &rho_left,
    const double &u_left,
    const double &p_left,
    const double &rho_right,
    const double &u_right,
    const double &p_right,
    const double &gamma,
    const double &x_interface
  );

  double value(
    const Point<dim>   &p,
    const unsigned int component) const override;

private:

  const double rho_left;
  const double u_left;
  const double p_left;
  const double rho_right;
  const double u_right;
  const double p_right;
  const double gamma;
  const double x_interface;
  const double g1;
  const double g2;
  const double g3;
  const double g4;
  const double g5;
  const double g6;
  const double g7;
  const double g8;
  const double c_left;
  const double c_right;
  double u_star;
  double p_star;

};

#include "EulerRiemannSolver.cc"
#endif
