/** \brief computes the exact solution for test problem 4 for a given x-point
 */
template <>
double ExactSolution4<1>::value (const Point<1>  &p,
                                 const unsigned int) const
{
  double return_value = p[0] * std::sin(p[0]) * std::exp(-p[0]);
  return return_value;
}
