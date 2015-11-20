/**
 * \file ShallowWaterLowOrderViscosity.h
 * \brief Provides the header for the ShallowWaterLowOrderViscosity class.
 */
#ifndef ShallowWaterLowOrderViscosity_h
#define ShallowWaterLowOrderViscosity_h

#include "include/viscosity/LowOrderViscosity.cc"

using namespace dealii;

/**
 * \class ShallowWaterLowOrderViscosity
 * \brief Class for low-order viscosity for the shallow water equations;
 *        this low-order viscosity allows the Froude number to be used.
 */
template <int dim>
class ShallowWaterLowOrderViscosity : public LowOrderViscosity<dim>
{
public:
  ShallowWaterLowOrderViscosity();
};

#include "src/viscosity/ShallowWaterLowOrderViscosity.cc"

#endif
