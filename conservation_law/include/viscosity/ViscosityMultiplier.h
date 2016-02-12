/**
 * \file ViscosityMultiplier.h
 * \brief Provides the header for the ViscosityMultiplier class.
 */
#ifndef ViscosityMultiplier_h
#define ViscosityMultiplier_h

#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

/**
 * \class ViscosityMultiplier
 * \brief Base class for computing multipliers for viscosity.
 *
 * This base class uses a default value of one for the multipliers;
 * derived classes override the `get_multiplier()` function.
 */
template <int dim>
class ViscosityMultiplier
{
public:
  ViscosityMultiplier();

  virtual double get_multiplier(const FEValues<dim> & fe_values,
                                const Vector<double> & solution) const;
};

#include "src/viscosity/ViscosityMultiplier.cc"

#endif
