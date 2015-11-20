/**
 * \file GraphTheoreticDiffusion.h
 * \brief Provides the header for the GraphTheoreticDiffusion class.
 */
#ifndef GraphTheoreticDiffusion_h
#define GraphTheoreticDiffusion_h

#include "include/diffusion/ArtificialDiffusion.h"
#include "include/viscosity/Viscosity.h"

using namespace dealii;

/**
 * \class GraphTheoreticDiffusion
 * \brief Class for graph-theoretic diffusion.
 */
template <int dim>
class GraphTheoreticDiffusion : public ArtificialDiffusion<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  GraphTheoreticDiffusion();

  void apply(std::shared_ptr<Viscosity<dim>> viscosity,
                     const Vector<double> & solution,
                     const Cell & cell,
                     const FEValues<dim> & fe_values,
                     Vector<double> & cell_residual) const override;
};

#include "src/diffusion/GraphTheoreticDiffusion.cc"

#endif
