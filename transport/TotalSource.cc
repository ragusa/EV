/** \brief computes the total source at a given point in space
 */
template<int dim>
double TotalSource<dim>::value(const Point<dim> &p) const
{
   double return_value = 0.0;
   switch (parameters.source_option) {
      case 1: {
         return_value = parameters.source_value / (4.0 * numbers::PI); // isotropic source term
         break;
      } case 2: {
         bool in_nonzero_region = true;
         for (unsigned int i = 0; i < dim; ++i)
            if (p[i] < 0.0) {
               in_nonzero_region = false;
               break;
            }
         if (in_nonzero_region)
            return_value = parameters.source_value;
         break;
      } case 3: {
         /* manufactured solution source
          * solution is: x*sin(x)*exp(-x) and sigma = 1
          */
         return_value = std::sin(p[0]) * std::exp(-p[0])
                         + p[0] * std::cos(p[0]) * std::exp(-p[0]);
         break;
      }
      default: {
         Assert(false,ExcNotImplemented());
         break;
      }
   }
   return return_value;
}

/** \brief computes the total source at a number of points in space
 */
template<int dim>
void TotalSource<dim>::value_list(const std::vector<Point<dim> > &points,
                                  std::vector<double> &values) const
{
   Assert(values.size() == points.size(),
         ExcDimensionMismatch (values.size(), points.size()));

   for (unsigned int i = 0; i < points.size(); ++i)
      values[i] = value(points[i]);
}
