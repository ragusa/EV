#ifndef ExactSolution1_cc
#define ExactSolution1_cc

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

using namespace dealii;

/** \brief defines the exact solution for test problem 1
 */
template <int dim>
class ExactSolution1 : public Function<dim> {
   public:
      ExactSolution1 () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int component = 0) const;
};

#include "ExactSolution1.cc"
#endif
