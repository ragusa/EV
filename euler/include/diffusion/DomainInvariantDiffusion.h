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
    const Triangulation<dim> & triangulation_,
    const QGauss<dim> & cell_quadrature_,
    const unsigned int & n_components,
    const unsigned int & n_dofs);

  void reinitialize();

  void compute_diffusion_matrix(const Vector<double> & solution,
                                const std::shared_ptr<Viscosity<dim>> viscosity,
                                SparseMatrix<double> & diffusion_matrix) override;

protected:
  void compute_gradients_and_normals();

  /** \brief Pointer to max wave speed object */
  const std::shared_ptr<MaxWaveSpeed<dim>> max_wave_speed;

  /** \brief Scalar 1st-order Lagrangian finite element */
  const FE_Q<dim> fe_scalar;

  /** \brief Degree of freedom handler for scalar case */
  DoFHandler<dim> dof_handler_scalar;

  /** \brief Pointer to cell quadrature */
  const QGauss<dim> * const cell_quadrature;

  /** \brief Number of solution components */
  const unsigned int n_components;

  /** \brief Number of degrees of freedom for the vector case */
  const unsigned int n_dofs;

  /** \brief Number of degrees of freedom for the scalar case */
  const unsigned int n_dofs_scalar;

  /** \brief Number of degrees of freedom per cell for the scalar case */
  const unsigned int dofs_per_cell_scalar;

  /** \brief Number of quadrature points per cell */
  const unsigned int n_q_points_cell;

  /** \brief Sparsity pattern for gradient matrices */
  SparsityPattern sparsity;

  /** \brief Array of matrices for each component of gradient tensor matrix */
  SparseMatrix<double> gradients[dim];

  /** \brief Sparse matrix of L2-norm of each C matrix entry */
  SparseMatrix<double> gradient_norms;
};

#include "src/diffusion/DomainInvariantDiffusion.cc"

#endif
