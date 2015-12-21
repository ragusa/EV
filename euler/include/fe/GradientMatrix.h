/**
 * \file GradientMatrix.h
 * \brief Provides the header for the GradientMatrix class.
 */
#ifndef GradientMatrix_h
#define GradientMatrix_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

/**
 * \class GradientMatrix
 * \brief Class for computing gradient matrix entries.
 */
template <int dim>
class GradientMatrix
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  GradientMatrix(const Triangulation<dim> & triangulation_,
                 const QGauss<dim> & cell_quadrature_,
                 const unsigned int & n_components);

  void reinitialize();

  Tensor<1, dim> get_normal(const unsigned int & i, const unsigned int & j) const;

  std::vector<Tensor<1, dim>> get_normals(const unsigned int & i) const;

  double get_gradient_norm(const unsigned int & i, const unsigned int & j) const;

  std::vector<double> get_gradient_norms(const unsigned int & i) const;

  std::vector<unsigned int> get_column_indices(const unsigned int & i) const;

protected:
  void compute_gradients_and_normals();

  /** \brief Scalar 1st-order Lagrangian finite element */
  const FE_Q<dim> fe_scalar;

  /** \brief Degree of freedom handler for scalar case */
  DoFHandler<dim> dof_handler_scalar;

  /** \brief Pointer to cell quadrature */
  const QGauss<dim> * const cell_quadrature;

  /** \brief Number of solution components */
  const unsigned int n_components;

  /** \brief Number of degrees of freedom for the vector case */
  unsigned int n_dofs;

  /** \brief Number of degrees of freedom for the scalar case */
  unsigned int n_dofs_scalar;

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

#include "src/fe/GradientMatrix.cc"

#endif
