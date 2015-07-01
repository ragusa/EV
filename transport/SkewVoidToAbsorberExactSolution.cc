/** \brief Constructor.
 */
template<int dim>
SkewVoidToAbsorberExactSolution<dim>::SkewVoidToAbsorberExactSolution() :
   Function<dim>(),
   transport_direction({1.0/sqrt(2.0),1.0/sqrt(3.0),1.0/sqrt(6.0)}),
   sigma(10.0),
   incoming(1.0)
{
}

/** \brief Computes value at a single point.
 */
template<int dim>
double SkewVoidToAbsorberExactSolution<dim>::value(
   const Point<dim>   &p,
   const unsigned int component) const
{
   // get time
   double t = this->get_time();

   // determine if x-Omega*t lies within domain; if not,
   // then solution hasn't reached x yet
   bool not_in_domain = false;
   for (unsigned int d = 0; d < dim; ++d)
      if (p[d] - transport_direction[d]*t < 0.0)
         not_in_domain = true;

   // compute exact solution
   double exact_solution;
   if (not_in_domain) {

      // determine if point is in absorber region
      bool in_absorber = true;
      for (unsigned int d = 0; d < dim; ++d)
         if (p[d] < 0.5)
            in_absorber = false;

      if (in_absorber) {

         // compute distance travelled in absorber
         double s = computeDistanceTravelledInAbsorber(p);

         // solution has been attenuated in absorber region
         exact_solution = incoming*std::exp(-sigma*s);

      } else {

         // solution has not been attenuated
         exact_solution = incoming;
      }

   } else {
      // solution has not reached x yet
      exact_solution = 0.0;
   }

   return exact_solution;
}

/** \brief Given that a point is in the absorber region, determines the distance
 *         travelled in absorber region.
 */
template<int dim>
double SkewVoidToAbsorberExactSolution<dim>::computeDistanceTravelledInAbsorber(
   const Point<dim> &p) const
{
   // components of distance travelled in absorber region
   Tensor<1,3> s;

   // index for direction that has the minimum distance ratio
   unsigned int d_min = 0;

   // initialize to x-direction
   s[0] = p[0] - 0.5;
   double r_min = s[0]/transport_direction[0];

   // loop over the number of dimensions and compare distance ratios
   // to determine which absorber region plane segment through which the
   // particle passes
   for (unsigned int d = 1; d < dim; ++d)
   {
      // compute distance ratio for this direction
      double r_d = (p[d] - 0.5)/transport_direction[d];

      // compare to the minimum distance comp
      if (r_d < r_min) { // this direction is new minimum distance ratio

         // set minimum distance ratio index to this direction
         d_min = d;

         // set minimum distance ratio to this one
         r_min = r_d;

         // compute distance for this direction
         s[d] = p[d] - 0.5;

         // recompute distance for all previous directions
         for (unsigned int dd = 0; dd <= d-1; ++dd)
            s[dd] = transport_direction[dd] / transport_direction[d] * s[d];         

      } else {
         // compute distance for this direction
         s[d] = transport_direction[d] / transport_direction[d_min] * s[d_min];
      }
   }

   // if number of dimensions for this problem is less than 3, compute
   // distance components for remaining dimensions up to 3
   for (unsigned int d = dim; d < 3; ++d)
      s[d] = transport_direction[d] / transport_direction[d_min] * s[d_min];

   // compute total distance
   return std::sqrt(std::pow(s[0],2) + std::pow(s[1],2) + std::pow(s[2],2));
}
