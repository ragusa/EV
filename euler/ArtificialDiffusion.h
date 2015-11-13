/**
 * \file ArtificialDiffusion.h
 * \brief Provides the header for the ArtificialDiffusion class.
 */
#ifndef ArtificialDiffusion_h
#define ArtificialDiffusion_h

using namespace dealii;

/**
 * \class ArtificialDiffusion
 * \brief Class for ...
 */
template <int dim>
class ArtificialDiffusion
{
public:
  ArtificialDiffusion();

  virtual void apply() const = 0;
};

#include "ArtificialDiffusion.cc"

#endif
