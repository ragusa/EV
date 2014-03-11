#ifndef ExactSolution4_cc
#define ExactSolution4_cc

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

using namespace dealii;

/** \brief defines the exact solution for test problem 4
 */
template <int dim>
class ExactSolution4 : public Function<dim> {
   public:
      ExactSolution4 () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int component = 0) const;
};

#include "ExactSolution4.cc"
#endif
