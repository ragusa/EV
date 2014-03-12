/** \brief computes the initial values
 */
template <>
double InitialValues<1>::value (const Point<1>  &p,
                                const unsigned int) const
{
  const double sigma  = 1.0;
  const double source = 1.0 / (4.0 * numbers::PI);

  double return_value;
  return_value = source*(1.0 - std::exp(-sigma*(p[0]+1.0)));

  return return_value;
}
