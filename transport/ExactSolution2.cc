/** \brief computes the exact solution for test problem 2 for a given x-point
 */
template <>
double ExactSolution2<1>::value (const Point<1>  &p,
                                 const unsigned int) const
{
  const double sigma = 100.0;
  const double incoming = 1.0;

  double return_value;
  if (p[0] < 0.0) return_value = incoming;
  else            return_value = incoming * std::exp(-sigma*p[0]);

  return return_value;
}
