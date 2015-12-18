/**
 * \file NoDiffusion.h
 * \brief Provides the header for the NoDiffusion class.
 */
#ifndef NoDiffusion_h
#define NoDiffusion_h

#include "include/diffusion/ArtificialDiffusion.h"

using namespace dealii;

/**
 * \class NoDiffusion
 * \brief Class for zero artificial diffusion. Does nothing.
 */
template <int dim>
class NoDiffusion : public ArtificialDiffusion<dim>
{
public:
  using Cell = typename ArtificialDiffusion<dim>::Cell;

  NoDiffusion();

  void compute_diffusion_matrix(const Vector<double> & solution,
                                const std::shared_ptr<Viscosity<dim>> viscosity,
                                SparseMatrix<double> & diffusion_matrix) override;
};

#include "src/diffusion/NoDiffusion.cc"

#endif
