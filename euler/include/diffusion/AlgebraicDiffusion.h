/**
 * \file AlgebraicDiffusion.h
 * \brief Provides the header for the AlgebraicDiffusion class.
 */
#ifndef AlgebraicDiffusion_h
#define AlgebraicDiffusion_h

#include "include/diffusion/ArtificialDiffusion.h"

using namespace dealii;

/**
 * \class AlgebraicDiffusion
 * \brief Class for artificial diffusion that is expressed using a diffusion
 *        matrix instead of artificial viscosities.
 */
template <int dim>
class AlgebraicDiffusion : public ArtificialDiffusion<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  AlgebraicDiffusion();

  void reinitialize(const SparsityPattern & sparsity,
                    const DoFHandler<dim> & dof_handler) override;

  void apply(std::shared_ptr<Viscosity<dim>> viscosity,
             const Vector<double> & solution,
             const Cell & cell,
             const FEValues<dim> & fe_values,
             Vector<double> & cell_residual) const override;

  void apply_algebraic_diffusion(const Vector<double> & solution,
                                 Vector<double> & ss_residual) override;

protected:
  /** \brief Temporary vector for matrix-vector product results */
  Vector<double> tmp_vector;
};

#include "src/diffusion/AlgebraicDiffusion.cc"

#endif
