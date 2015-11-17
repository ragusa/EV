/**
 * \file NoDiffusion.h
 * \brief Provides the header for the NoDiffusion class.
 */
#ifndef NoDiffusion_h
#define NoDiffusion_h

#include "ArtificialDiffusion.h"

using namespace dealii;

/**
 * \class NoDiffusion
 * \brief Class for having no artificial diffusion. Does nothing.
 */
template <int dim>
class NoDiffusion : public ArtificialDiffusion<dim>
{
public:
  using Cell = typename ArtificialDiffusion<dim>::Cell;

  NoDiffusion();

  void apply(std::shared_ptr<Viscosity<dim>> viscosity,
             const Vector<double> & solution,
             const Cell & cell,
             const FEValues<dim> & fe_values,
             Vector<double> & cell_residual) const override;
};

#include "NoDiffusion.cc"

#endif
