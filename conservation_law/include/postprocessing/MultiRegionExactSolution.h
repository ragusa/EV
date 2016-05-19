#ifndef MultiRegionExactSolution_cc
#define MultiRegionExactSolution_cc

#include <deal.II/base/function.h>

using namespace dealii;

/**
 * \brief Exact solution function class for the multi-region unit hypercube
 *        test problem.
 */
template <int dim>
class MultiRegionExactSolution : public Function<dim>
{
public:
  MultiRegionExactSolution(const std::vector<double> & interface_positions,
                           const std::vector<double> & source,
                           const std::vector<double> & sigma,
                           const Tensor<1, dim> & direction,
                           const double & incoming);

  double value(const Point<dim> & p,
               const unsigned int component = 0) const override;

protected:
  /** number of regions */
  const unsigned int Nr;

  /** "radii" between each region, including boundary radii */
  std::vector<double> ri;

  /** cross sections for each region */
  const std::vector<double> sigma;

  /** source for each region */
  const std::vector<double> q;

  /** transport speed */
  const double c;

  /** transport direction */
  const Tensor<1, dim> direction;

  /** incoming flux value */
  const double incoming;

private:
  double computeDistanceTravelledInRegion(const Point<dim> & p,
                                          const double & r) const;

  virtual void adjust_distances(std::vector<double> & s) const = 0;

  virtual double compute_reference_solution(const Point<dim> & p) const = 0;
};

#include "src/postprocessing/MultiRegionExactSolution.cc"
#endif
