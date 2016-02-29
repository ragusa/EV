/**
 * \file DomainInvariantDiffusion.h
 * \brief Provides the header for the DomainInvariantDiffusion class.
 */
#ifndef DomainInvariantDiffusion_h
#define DomainInvariantDiffusion_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "include/diffusion/ArtificialDiffusion.h"
#include "include/fe/GradientMatrix.h"
#include "include/viscosity/MaxWaveSpeed.h"

using namespace dealii;

/**
 * \class DomainInvariantDiffusion
 * \brief Class for a domain-invariant diffusion matrix.
 *
 * This matrix is computed as
 * \f[
 *   D_{i,j} = \left\{\begin{array}{c l}
 *     -\max\left(
 *     \lambda_{max}(\mathbf{n}_{i,j},\mathbf{u}_i^n,\mathbf{u}_j^n)
 *     \|\mathbf{c}_{i,j}\|_{\ell^2},
 *     \lambda_{max}(\mathbf{n}_{j,i},\mathbf{u}_j^n,\mathbf{u}_i^n)
 *     \|\mathbf{c}_{j,i}\|_{\ell^2}
 *     \right) & j\ne i \\
 *     -\sum\limits_{k\ne i}D_{i,k} & j = i
 *     \end{array}\right. \,.
 * \f]
 */
template <int dim>
class DomainInvariantDiffusion : public ArtificialDiffusion<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename ArtificialDiffusion<dim>::Cell;

  DomainInvariantDiffusion(
    const std::shared_ptr<MaxWaveSpeed<dim>> & max_wave_speed_,
    const std::shared_ptr<GradientMatrix<dim>> & gradient_matrix,
    const unsigned int & n_components,
    const unsigned int & n_dofs);

  void reinitialize();

  void compute_diffusion_matrix(const Vector<double> & solution,
                                const std::shared_ptr<Viscosity<dim>> viscosity,
                                SparseMatrix<double> & diffusion_matrix) override;

protected:
  /** \brief Pointer to max wave speed object */
  const std::shared_ptr<MaxWaveSpeed<dim>> max_wave_speed;

  /** \brief Pointer to gradient matrix object */
  const std::shared_ptr<GradientMatrix<dim>> gradient_matrix;

  /** \brief Number of solution components */
  const unsigned int n_components;

  /** \brief Number of degrees of freedom for the scalar case */
  const unsigned int n_dofs_scalar;
};

#include "src/diffusion/DomainInvariantDiffusion.cc"

#endif
