#ifndef TransientMultiRegionExactSolution_cc
#define TransientMultiRegionExactSolution_cc

#include "include/postprocessing/MultiRegionExactSolution.h"

using namespace dealii;

/**
 * \brief Exact solution function class for the multi-region unit hypercube
 *        transient test problem.
 */
template <int dim>
class TransientMultiRegionExactSolution : public MultiRegionExactSolution<dim>
{
public:
  TransientMultiRegionExactSolution(
    const std::vector<double> & interface_positions,
    const std::vector<double> & source,
    const std::vector<double> & sigma,
    const Tensor<1, dim> & direction,
    const double & incoming);

private:
  void adjust_distances(std::vector<double> & s) const override;

  double compute_reference_solution(const Point<dim> & p) const override;
};

#include "src/postprocessing/TransientMultiRegionExactSolution.cc"
#endif
