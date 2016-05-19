#ifndef SteadyStateMultiRegionExactSolution_cc
#define SteadyStateMultiRegionExactSolution_cc

#include "include/postprocessing/MultiRegionExactSolution.h"

using namespace dealii;

/**
 * \brief Exact solution function class for the multi-region unit hypercube
 *        steady-state test problem.
 */
template <int dim>
class SteadyStateMultiRegionExactSolution : public MultiRegionExactSolution<dim>
{
public:
  SteadyStateMultiRegionExactSolution(
    const std::vector<double> & interface_positions,
    const std::vector<double> & source,
    const std::vector<double> & sigma,
    const Tensor<1, dim> & direction,
    const double & incoming);

private:
  void adjust_distances(std::vector<double> & s) const override;

  double compute_reference_solution(const Point<dim> & p) const override;
};

#include "src/postprocessing/SteadyStateMultiRegionExactSolution.cc"
#endif
