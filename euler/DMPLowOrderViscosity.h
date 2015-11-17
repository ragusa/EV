/**
 * \file DMPLowOrderViscosity.h
 * \brief Provides the header for the DMPLowOrderViscosity class.
 */
#ifndef DMPLowOrderViscosity_h
#define DMPLowOrderViscosity_h

using namespace dealii;

/**
 * \class DMPLowOrderViscosity
 * \brief Class for low-order viscosity that satisfies a discrete maximum
 *        principle. This is only valid for scalar conservation laws.
 */
template <int dim>
class DMPLowOrderViscosity : public Viscosity<dim>
{
public:
  DMPLowOrderViscosity();
};

#include "DMPLowOrderViscosity.cc"

#endif
