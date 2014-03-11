#ifndef TotalSource_cc
#define TotalSource_cc

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include "TransportParameters.h"

/** \brief defines the total source function
 */
template<int dim>
class TotalSource: public Function<dim>
{
   public:
      TotalSource(const TransportParameters &parameters) :
         Function<dim>(),
         parameters(parameters)
      {}

      double value(const Point<dim> &p) const;

      void value_list(const std::vector<Point<dim> > &points,
                      std::vector<double> &values) const;
   private:
      const TransportParameters &parameters;
};

#include "TotalSource.cc"
#endif
