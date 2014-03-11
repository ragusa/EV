#ifndef ExactSolution2_cc
#define ExactSolution2_cc

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

using namespace dealii;

/** \brief defines the exact solution for test problem 2
 */
template <int dim>
class ExactSolution2 : public Function<dim> {
   public:
      ExactSolution2 () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int component = 0) const;
};

#include "ExactSolution2.cc"
#endif
