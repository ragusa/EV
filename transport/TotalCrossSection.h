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
      TotalCrossSection(const unsigned int &cross_section_option, const double &cross_section_value) :
         Function<dim>(),
         cross_section_option(cross_section_option),
         cross_section_value(cross_section_value)
      {}

      double value(const Point<dim> &p) const;

      void value_list(const std::vector<Point<dim> > &points,
                      std::vector<double> &values) const;
   private:
      const unsigned int cross_section_option;
      const double cross_section_value;
};

#include "TotalCrossSection.cc"
#endif
