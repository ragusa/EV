/**
 * \file ShallowWaterEntropy.h
 * \brief Provides the header for the ShallowWaterEntropy class.
 */
#ifndef ShallowWaterEntropy_h
#define ShallowWaterEntropy_h

#include "ShallowWaterEntropyFluxFEValuesCell.h"
#include "ShallowWaterEntropyFluxFEValuesFace.h"

using namespace dealii;

/**
 * \class ShallowWaterEntropy
 * \brief Class for entropy for the shallow water equations.
 */
template <int dim>
class ShallowWaterEntropy : public Entropy<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Entropy<dim>::Cell;

  ShallowWaterEntropy();

  std::vector<double> compute_entropy() const override;

  std::vector<double> compute_divergence_entropy_flux(
    const Cell & cell) const override;

  void reinitialize_group_fe_values(const Vector<double> & solution) override;

private:
  const FEValuesExtractors::Scalar height_extractor;

  const FEValuesExtractors::Vector momentum_extractor;

  const double gravity;

  ShallowWaterEntropyFluxFEValuesCell<dim> entropy_flux_fe_values_cell;

  ShallowWaterEntropyFluxFEValuesFace<dim> entropy_flux_fe_values_face;
};

#include "ShallowWaterEntropy.cc"

#endif
