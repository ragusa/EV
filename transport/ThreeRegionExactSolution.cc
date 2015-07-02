/** Constructor.
 */
template<int dim>
ThreeRegionExactSolution<dim>::ThreeRegionExactSolution() :
   Function<dim>(),
   Nr(3),
   xi({0.0,0.3,0.6,1.0}),
   sigma({1.0,40.0,20.0}),
   q({1.0,5.0,20.0}),
   c(1.0),
   incoming(1.0)
{
}

/** Computes exact solution at a single point.
 */
template<int dim>
double ThreeRegionExactSolution<dim>::value(
   const Point<dim>   &p,
   const unsigned int component) const
{
   // get time
   double t = this->get_time();

   // get position
   double x = p[0];

   // get reference position x-c*t
   double xref = x - c*t;

   // compute distances traversed in each region
   std::vector<double> s(Nr);
   for (unsigned int i = 0; i < Nr; ++i) {
      double splus  = std::min(x,xi[i+1]);
      double sminus = std::max(xref,xi[i]);
      s[i] = std::max(splus-sminus,0.0);
   }

   // compute solution component due to boundary flux
   // initialize to reference solution
   double ub;
   if (xref < xi[0])
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

