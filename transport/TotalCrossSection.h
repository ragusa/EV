#ifndef TotalCrossSection_cc
#define TotalCrossSection_cc

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include "TransportParameters.h"

/** \brief defines the total cross section function
 */
template<int dim>
class TotalCrossSection: public Function<dim>
{
   public:
      TotalCrossSection(const TransportParameters &parameters) :
         Function<dim>(),
         parameters(parameters)
      {}

      double value(const Point<dim> &p) const;

      void value_list(const std::vector<Point<dim> > &points,
                      std::vector<double> &values) const;
   private:
      const TransportParameters parameters;
};

#include "TotalCrossSection.cc"
#endif
