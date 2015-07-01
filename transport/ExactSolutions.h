#ifndef ExactSolutions_cc
#define ExactSolutions_cc

#include <deal.II/base/function.h>

using namespace dealii;

enum ExactSolutionOption {none, parser, skew_void_to_absorber};

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

#include "SkewVoidToAbsorberExactSolution.cc"
#endif
