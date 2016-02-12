/**
 * \file ShallowWaterViscosityMultiplier.h
 * \brief Provides the header for the ShallowWaterViscosityMultiplier class.
 */
#ifndef ShallowWaterViscosityMultiplier_h
#define ShallowWaterViscosityMultiplier_h

#include <deal.II/fe/fe_values_extractors.h>

#include "include/viscosity/ViscosityMultiplier.h"

using namespace dealii;

/**
 * \class ShallowWaterViscosityMultiplier
 * \brief Class for computing Froude number to be used as multiplier for
 *        viscosity.
 */
template <int dim>
class ShallowWaterViscosityMultiplier : public ViscosityMultiplier<dim>
{
public:
  ShallowWaterViscosityMultiplier(const double & gravity);

  double get_multiplier(const FEValues<dim> & fe_values,
                        const Vector<double> & solution) const override;

private:
  /** \brief Acceleration due to gravity */
  const double gravity;

  /** \brief FE values extractor for height */
  const FEValuesExtractors::Scalar height_extractor;

  /** \brief FE values extractor for momentum */
  const FEValuesExtractors::Vector momentum_extractor;
};

#include "src/viscosity/ShallowWaterViscosityMultiplier.cc"

#endif
