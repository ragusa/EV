/**
 * \file ShallowWaterPostProcessor.h
 * \brief Contains the header for the ShallowWaterPostProcessor class.
 */

#ifndef ShallowWaterPostProcessor_cc
#define ShallowWaterPostProcessor_cc

#include <deal.II/numerics/data_postprocessor.h>

using namespace dealii;

/**
 * \brief Allows for the output of derived quantities related to the shallow
 *        water equations.
 */
template <int dim>
class ShallowWaterPostProcessor : public DataPostprocessor<dim>
{
public:
  ShallowWaterPostProcessor(
    const std::shared_ptr<Function<dim>> & bathymetry_function);

  std::vector<std::string> get_names() const override;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation() const override;

  UpdateFlags get_needed_update_flags() const override;

  void compute_derived_quantities_vector(
    const std::vector<Vector<double>> & uh,
    const std::vector<std::vector<Tensor<1, dim>>> & duh,
    const std::vector<std::vector<Tensor<2, dim>>> & dduh,
    const std::vector<Point<dim>> & normals,
    const std::vector<Point<dim>> & evaluation_points,
    std::vector<Vector<double>> & computed_quantities) const override;

private:
  std::shared_ptr<Function<dim>> bathymetry_function;
};

#include "src/postprocessing/ShallowWaterPostProcessor.cc"

#endif
