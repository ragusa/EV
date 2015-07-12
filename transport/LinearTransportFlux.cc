template<int dim>
LinearTransportFlux<dim>::LinearTransportFlux(const Tensor<1,dim> & velocity) :
    ConservationLawFlux<dim>(false),
    velocity(velocity)
{
}

template<int dim>
std::vector<Tensor<1,dim> > LinearTransportFlux<dim>::evaluate(
  const std::vector<double> & U,
  const std::vector<double> & x) const
{
  // get number of values
  const unsigned int n = U.size();

  // evaluate flux
  std::vector<Tensor<1,dim> > flux_values(n);
  for (unsigned int i = 0; i < n; ++i)
    flux_values[i] = U[i]*velocity;
  
  return flux_values;
}

template<int dim>
std::vector<Tensor<1,dim> > LinearTransportFlux<dim>::evaluate_derivative(
  const std::vector<double> & U,
  const std::vector<double> & x) const
{
  // get number of values
  const unsigned int n = U.size();

  // evaluate flux derivative
  std::vector<Tensor<1,dim> > flux_derivative_values(n,velocity);
  
  return flux_derivative_values;
}
