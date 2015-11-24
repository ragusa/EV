/**
 * \file DomainInvariantViscosity.h
 * \brief Provides the header for the DomainInvariantViscosity class.
 */
#ifndef DomainInvariantViscosity_h
#define DomainInvariantViscosity_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "include/viscosity/Viscosity.h"

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

  /** \brief Alias for line iterator */
  using Line = typename DoFHandler<dim>::line_iterator;

  DomainInvariantViscosity(const Triangulation<dim> & triangulation);

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n);

protected:
  virtual double compute_max_wave_speed(
    const Tensor<1, dim> normal,
    const std::vector<double> solution_i,
    const std::vector<double> solution_j) const = 0;

  void compute_graph_theoretic_sums();

  void compute_gradients_and_normals();

  /** \brief Scalar 1st-order Lagrangian finite element */
  const FE_Q<dim> fe;

  /** \brief Degree of freedom handler for a scalar */
  DoFHandler<dim> dof_handler;

  /** \brief Pointer to triangulation */
  const Triangulation<dim> * const triangulation;

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

  /** \brief Map of line iterator to max viscosity along line */
  std::map<Line, double> max_viscosity;

  /** \brief Number of lines in triangulation */
  unsigned int n_lines;

  /** \brief Vector of gradient tensors from 1st node on each line */
  std::vector<Tensor<1, dim>> c_ij;

  /** \brief Vector of gradient tensors from 2nd node on each line */
  std::vector<Tensor<1, dim>> c_ji;

  /** \brief Vector of norm of gradient tensors from 1st node on each line */
  std::vector<double> c_ij_norm;

  /** \brief Vector of norm of gradient tensors from 2nd node on each line */
  std::vector<double> c_ji_norm;

  /** \brief Vector of normal vectors from 1st node on each line */
  std::vector<Tensor<1, dim>> normal_ij;

  /** \brief Vector of normal vectors from 2nd node on each line */
  std::vector<Tensor<1, dim>> normal_ji;
};

#include "src/viscosity/DomainInvariantViscosity.cc"

#endif
