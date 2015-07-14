#ifndef ExactSolutions_cc
#define ExactSolutions_cc

#include <deal.II/base/function.h>

using namespace dealii;

/** Enumeration for type of exact solution */
enum ExactSolutionOption
{
  none,
  parser,
  multi_region
};

/** Exact solution function class for the multi-region unit hypercube
 *  test problem.
 */
template<int dim>
class MultiRegionExactSolution : public Function<dim> {
   public:
      MultiRegionExactSolution(
        const std::vector<double> &interface_positions,
        const std::vector<double> &source,
        const std::vector<double> &sigma,
        const std::vector<double> &direction,
        const double              &incoming
      );

      double value(const Point<dim>   &p,
                   const unsigned int component = 0) const override;

   private:

      double computeDistanceTravelledInRegion(
         const Point<dim> &p,
         const double     &r) const;

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
      const std::vector<double> direction;

      /** incoming flux value */
      const double incoming;
};

#include "MultiRegionExactSolution.cc"
#endif
