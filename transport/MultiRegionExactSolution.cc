/** Constructor.
 */
template<int dim>
MultiRegionExactSolution<dim>::MultiRegionExactSolution(
  const std::vector<double> &interface_positions,
  const std::vector<double> &source,
  const std::vector<double> &sigma,
  const std::vector<double> &direction,
  const double              &incoming
) :
    Function<dim>(),
    Nr(sigma.size()),
    sigma(sigma),
    q(source),
    c(1.0),
    direction(direction),
    incoming(incoming)
{
   ri.resize(Nr+1);
   ri[0] = 0.0;
   for (unsigned int i = 1; i <= Nr-1; ++i)
      ri[i] = interface_positions[i-1];
   ri[Nr] = 1.0; 
}

/** Computes exact solution at a single point.
 */
template<int dim>
double MultiRegionExactSolution<dim>::value(
   const Point<dim>   &p,
   const unsigned int component) const
{
   // get time
   double t = this->get_time();

   // copy point
   Point<dim> x(p);

   // compute distances traversed in each region
   std::vector<double> s(Nr);
   for (int i = Nr-1; i >= 0; --i)
   {
      // determine if point is in region
      bool entered_region = true;
      for (unsigned int d = 0; d < dim; ++d)
         if (x[d] < ri[i])
            entered_region = false;

      // if in the region, compute the distance in the region
      if (entered_region) {

         // compute distance in the region
         s[i] = computeDistanceTravelledInRegion(x,ri[i]);

         // shift point to the entrance point for computation of
         // previous region distance
         for (unsigned int d = 0; d < dim; ++d)
            x[d] -= s[i]*direction[d];
         
      } else {
         s[i] = 0.0;
      }
   }

   // the above computations assume that x-c*Omega*t is outside the
   // domain, so distances need to be adjusted if this is not true
   // first, compute sum of distances
   double s_sum = 0.0;
   for (unsigned int i = 0; i < Nr; ++i)
      s_sum += s[i];
   // compute excess distance
   double excess_distance = s_sum - c*t;
   // if there is an excess, adjust distances
   if (excess_distance > 0.0)
   {
      for (unsigned int i = 0; i < Nr; ++i)
      {
         if (excess_distance > s[i])
         {
            excess_distance -= s[i];
            s[i] = 0.0;
         } else {
            s[i] -= excess_distance;
            excess_distance = 0.0;
            break;
         }
      }
   }

   // determine if x-c*Omega*t lies within domain; if not,
   // then solution hasn't reached x yet
   bool not_in_domain = false;
   for (unsigned int d = 0; d < dim; ++d)
      if (p[d] - direction[d]*c*t < 0.0)
         not_in_domain = true;

   // compute solution component due to boundary flux
   // initialize to reference solution
   double ub;
   if (not_in_domain)
      ub = incoming;
   else
      ub = 0.0;

   // compute the number of mean free paths from xref to x
   double ub_mfp = 0.0;
   for (unsigned int i = 0; i < Nr; ++i)
      ub_mfp += sigma[i]*s[i];
   // attenuate reference solution
   ub *= exp(-ub_mfp);

   // compute solution component due to sources
   double uq = 0.0;
   for (unsigned int i = 0; i < Nr; ++i) {
      // compute base source contribution from region
      double uqi;
      if (sigma[i] == 0.0)
         uqi = q[i]*s[i];
      else
         uqi = q[i]/sigma[i]*(1.0 - exp(-sigma[i]*s[i]));

      // compute the number of mean free paths in subsequent regions
      double uq_mfp = 0.0;
      for (unsigned int j = i+1; j < Nr; ++j)
         uq_mfp += sigma[j]*s[j];

      // add attenuated source to total
      uq += uqi*exp(-uq_mfp);
   }

   // return total solution
   return ub + uq;
}

/** Given that a point is in a region, determines the distance
 *  travelled in that region.
 */
template<int dim>
double MultiRegionExactSolution<dim>::computeDistanceTravelledInRegion(
  const Point<dim> &p,
  const double     &r
) const
{
   // components of distance travelled in region
   std::vector<double> s(3,0.0);

   // index for direction that has the minimum distance ratio
   unsigned int d_min = 0;

   // initialize to x-direction
   s[0] = p[0] - r;
   double r_min = s[0]/direction[0];

   // loop over the number of dimensions and compare distance ratios
   // to determine which absorber region plane segment through which the
   // particle passes
   for (unsigned int d = 1; d < dim; ++d)
   {
      // compute distance ratio for this direction
      double r_d = (p[d] - r)/direction[d];

      // compare to the minimum distance comp
      if (r_d < r_min) { // this direction is new minimum distance ratio

         // set minimum distance ratio index to this direction
         d_min = d;

         // set minimum distance ratio to this one
         r_min = r_d;

         // compute distance for this direction
         s[d] = p[d] - r;

         // recompute distance for all previous directions
         for (unsigned int dd = 0; dd <= d-1; ++dd)
            s[dd] = direction[dd] / direction[d] * s[d];

      } else {
         // compute distance for this direction
         s[d] = direction[d] / direction[d_min] * s[d_min];
      }
   }

   // if number of dimensions for this problem is less than 3, compute
   // distance components for remaining dimensions up to 3
   for (unsigned int d = dim; d < 3; ++d)
      s[d] = direction[d] / direction[d_min] * s[d_min];

   // compute total distance
   return std::sqrt(std::pow(s[0],2) + std::pow(s[1],2) + std::pow(s[2],2));
}
