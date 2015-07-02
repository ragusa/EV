#ifndef ExactSolutions_cc
#define ExactSolutions_cc

#include <deal.II/base/function.h>

using namespace dealii;

enum ExactSolutionOption {none, parser, skew_void_to_absorber, three_region};

/** Exact solution function class for the skew void-to-absorber test problem.
 */
template<int dim>
class SkewVoidToAbsorberExactSolution : public Function<dim> {
   public:
      SkewVoidToAbsorberExactSolution();

      double value(const Point<dim>   &p,
                   const unsigned int component = 0) const override;

   private:
      double computeDistanceTravelledInAbsorber(
         const Point<dim> &p) const;

      const std::vector<double> transport_direction;
      const double sigma;
      const double incoming;
};

/** Exact solution function class for the 3-region test problem.
 */
template<int dim>
class ThreeRegionExactSolution : public Function<dim> {
   public:
      ThreeRegionExactSolution();

      double value(const Point<dim>   &p,
                   const unsigned int component = 0) const override;

   private:

      /** number of regions */
      const unsigned int Nr;

      /** x-points between each region, including left and right endpoints */
      const std::vector<double> xi;

      /** cross sections for each region */
      const std::vector<double> sigma;

      /** source for each region */
      const std::vector<double> q;

      /** transport speed */
      const double c;

      /** incoming flux value */
      const double incoming;
};

#include "SkewVoidToAbsorberExactSolution.cc"
#include "ThreeRegionExactSolution.cc"
#endif
