/**
 * \file DomainInvariantViscosity.h
 * \brief Provides the header for the DomainInvariantViscosity class.
 */
#ifndef DomainInvariantViscosity_h
#define DomainInvariantViscosity_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "include/viscosity/MaxWaveSpeed.h"
#include "include/viscosity/Viscosity.h"
#include "include/viscosity/ViscosityMultiplier.h"

using namespace dealii;

/**
 * \class DomainInvariantViscosity
 * \brief Class for domain-invariant low-order viscosity.
 */
template <int dim>
class DomainInvariantViscosity : public Viscosity<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Viscosity<dim>::Cell;

  DomainInvariantViscosity(
    const std::shared_ptr<MaxWaveSpeed<dim>> & max_wave_speed_,
    const FESystem<dim> & fe,
    const DoFHandler<dim> & dof_handler_,
    const Triangulation<dim> & triangulation_,
    const QGauss<dim> & cell_quadrature_,
    const unsigned int & n_components,
    const std::shared_ptr<ViscosityMultiplier<dim>> & viscosity_multiplier);

  void reinitialize() override;

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n) override;

protected:
  void compute_graph_theoretic_sums();

  void compute_gradients_and_normals();

  /** \brief Pointer to max wave speed object */
  const std::shared_ptr<MaxWaveSpeed<dim>> max_wave_speed;

  /** \brief Scalar 1st-order Lagrangian finite element */
  const FE_Q<dim> fe_scalar;

  /** \brief Pointer to finite element system */
  const FESystem<dim> * fe;

  /** \brief Degree of freedom handler for scalar case */
  DoFHandler<dim> dof_handler_scalar;

  /** \brief Degree of freedom handler for scalar case */
  const DoFHandler<dim> * const dof_handler;

  /** \brief Pointer to cell quadrature */
  const QGauss<dim> * const cell_quadrature;

  /** \brief Number of solution components */
  const unsigned int n_components;

  /** \brief Number of degrees of freedom for the scalar case */
  const unsigned int n_dofs_scalar;

  /** \brief Number of degrees of freedom per cell for the scalar case */
  const unsigned int dofs_per_cell_scalar;

  /** \brief Number of quadrature points per cell */
  const unsigned int n_q_points_cell;

  /** \brief Sparsity pattern for viscous bilinear form sum matrix */
  SparsityPattern sparsity;

  /**
   * \brief Matrix storing viscous bilinear form sums.
   *
   * Each element of this matrix is
   * \f[
   *   B_{i,j} = \sum_{K:D_K\subset S_{i,j}}b_K(\varphi_i,\varphi_j)
   * \f]
   */
  SparseMatrix<double> viscous_sums;

  /** \brief Array of matrices for each component of gradient tensor matrix */
  SparseMatrix<double> gradients[dim];

  /** \brief Sparse matrix of L2-norm of each C matrix entry */
  SparseMatrix<double> gradient_norms;

  /** \brief Pointer to viscosity multiplier */
  std::shared_ptr<ViscosityMultiplier<dim>> viscosity_multiplier;
};

#include "src/viscosity/DomainInvariantViscosity.cc"

#endif
