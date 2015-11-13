/**
 * \file ShallowWaterEntropy.h
 * \brief Provides the header for the ShallowWaterEntropy class.
 */
#ifndef ShallowWaterEntropy_h
#define ShallowWaterEntropy_h

using namespace dealii;

/**
 * \class ShallowWaterEntropy
 * \brief Class for entropy for the shallow water equations.
 */
template <int dim>
class ShallowWaterEntropy : public Entropy<dim>
{
public:
  ShallowWaterEntropy();
};

#include "ShallowWaterEntropy.cc"

#endif
