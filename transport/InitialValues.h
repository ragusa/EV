#ifndef InitialValues_cc
#define InitialValues_cc

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

using namespace dealii;

/** \brief defines the initial values
 */
template <int dim>
class InitialValues : public Function<dim> {
   public:
      InitialValues () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int component = 0) const;
};

#include "InitialValues.cc"
#endif
