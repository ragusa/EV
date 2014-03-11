/** \brief computes the total cross section at a given point in space
 */
template<int dim>
double TotalCrossSection<dim>::value(const unsigned int group,
                                     const unsigned int direction,
                                     const Point<dim> &p) const
{
   Assert(group < 2, ExcNotImplemented());
   double return_value = 0.0;
   switch (parameters.total_cross_section_option) {
      case 1: {
         if (group == 0)
            return_value = parameters.total_cross_section_value;
         else if (group == 1)
            return_value = 0.0;
         else {
            Assert(false,ExcNotImplemented());
         }
         break;
      }
      case 2: {
         if (group == 0) {
            bool in_nonzero_region = true;
            for (unsigned int i = 0; i < dim; ++i)
               if (p[i] < 0.0) {
                  in_nonzero_region = false;
                  break;
               }
            if (in_nonzero_region)
               return_value = parameters.total_cross_section_value;
         } else if (group == 1)
            return_value = 0.0;
         else {
            Assert(false,ExcNotImplemented());
         }
         break;
      }
      default: {
         Assert(false,ExcNotImplemented());
         break;
      }
   }

   return return_value;
}

/** \brief computes the total cross section at a number of points in space
 */
template<int dim>
void TotalCrossSection<dim>::value_list(const unsigned int group,
                                        const unsigned int direction,
                                        const std::vector<Point<dim> > &points,
                                        std::vector<double> &values) const
{
   Assert(values.size() == points.size(),
         ExcDimensionMismatch (values.size(), points.size()));

   for (unsigned int i = 0; i < points.size(); ++i)
      values[i] = TotalCrossSection<dim>::value(group, direction, points[i]);
}
