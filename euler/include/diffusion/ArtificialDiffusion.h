/**
 * \file ArtificialDiffusion.h
 * \brief Provides the header for the ArtificialDiffusion class.
 */
#ifndef ArtificialDiffusion_h
#define ArtificialDiffusion_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>

#include "include/viscosity/Viscosity.h"

using namespace dealii;

/**
 * \class ArtificialDiffusion
 * \brief Class for artificial diffusion
 */
template <int dim>
class ArtificialDiffusion
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  ArtificialDiffusion();

  virtual void reinitialize(const SparsityPattern & sparsity,
                            const DoFHandler<dim> & dof_handler);

  virtual void update(const Vector<double> & solution);

  virtual void apply(std::shared_ptr<Viscosity<dim>> viscosity,
                     const Vector<double> & solution,
                     const Cell & cell,
                     const FEValues<dim> & fe_values,
                     Vector<double> & cell_residual) const = 0;

  virtual void apply_algebraic_diffusion(const Vector<double> & solution,
                                         Vector<double> & ss_residual);

protected:
  /** \brief Diffusion matrix */
  SparseMatrix<double> diffusion_matrix;
};

#include "src/diffusion/ArtificialDiffusion.cc"

#endif
