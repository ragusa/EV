/** \brief computes the total cross section at a given point in space
 */
template<int dim>
double TotalCrossSection<dim>::value(const Point<dim> &p) const
{
   double return_value = 0.0;
   switch (cross_section_option) {
      case 1: {
         return_value = cross_section_value;
         break;
      }
      case 2: {
         bool in_nonzero_region = true;
         for (unsigned int i = 0; i < dim; ++i)
            if (p[i] < 0.0) {
               in_nonzero_region = false;
               break;
            }
         if (in_nonzero_region)
            return_value = cross_section_value;
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
void TotalCrossSection<dim>::value_list(const std::vector<Point<dim> > &points,
                                        std::vector<double> &values) const
{
   Assert(values.size() == points.size(),
         ExcDimensionMismatch (values.size(), points.size()));

   for (unsigned int i = 0; i < points.size(); ++i)
      values[i] = value(points[i]);
}
