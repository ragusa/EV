/**
 * \file Entropy.h
 * \brief Provides the header for the Entropy class.
 */
#ifndef Entropy_h
#define Entropy_h

using namespace dealii;

/**
 * \class Entropy
 * \brief Class for entropy.
 */
template <int dim>
class Entropy
{
public:
  Entropy();

private:
  virtual void compute_entropy() const = 0;
};

#include "Entropy.cc"

#endif
