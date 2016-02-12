/**
 * \file ShallowWaterPostProcessor.cc
 * \brief Provides the function definitions for the ShallowWaterPostProcessor class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ShallowWaterPostProcessor<dim>::ShallowWaterPostProcessor(
  const std::shared_ptr<Function<dim>> & bathymetry_function_)
  : bathymetry_function(bathymetry_function_)
{
}

/**
 * \brief Gets the names of the quantities included in the post-processor.
 */
template <int dim>
std::vector<std::string> ShallowWaterPostProcessor<dim>::get_names() const
{
  std::vector<std::string> names(2);
  names[0] = "bathymetry";
  names[1] = "waterlevel";
  return names;
}

/**
 * \brief Gets the interpretations of the quantities included in the
 *        post-processor.
 */
template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
  ShallowWaterPostProcessor<dim>::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretations(2, DataComponentInterpretation::component_is_scalar);
  return interpretations;
}

/**
 * \brief Gets the necessary update flags for evaluating post-processor
 *        quantities.
 */
template <int dim>
UpdateFlags ShallowWaterPostProcessor<dim>::get_needed_update_flags() const
{
  return update_values | update_q_points;
}

/**
 * \brief Computes the derived quantities at quadrature points.
 */
template <int dim>
void ShallowWaterPostProcessor<dim>::compute_derived_quantities_vector(
  const std::vector<Vector<double>> & uh,
  const std::vector<std::vector<Tensor<1, dim>>> &,
  const std::vector<std::vector<Tensor<2, dim>>> &,
  const std::vector<Point<dim>> &,
  const std::vector<Point<dim>> & evaluation_points,
  std::vector<Vector<double>> & computed_quantities) const
{
  const unsigned int n_quadrature_points = evaluation_points.size();
  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    const double b = bathymetry_function->value(evaluation_points[q]);
    computed_quantities[q](0) = b;            // bathymetry
    computed_quantities[q](1) = b + uh[q](0); // water level; b + h
  }
}
