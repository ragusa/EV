/**
 * \file ScalarEntropy.h
 * \brief Provides the header for the ScalarEntropy class.
 */
#ifndef ScalarEntropy_h
#define ScalarEntropy_h

#include "include/entropy/Entropy.h"

using namespace dealii;

/**
 * \class ScalarEntropy
 * \brief Class for entropy for a scalar conservation law equation.
 */
template <int dim>
class ScalarEntropy : public Entropy<dim>
{
public:
  ScalarEntropy();
};

#include "src/entropy/ScalarEntropy.cc"

#endif
