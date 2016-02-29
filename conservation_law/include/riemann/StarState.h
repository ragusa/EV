/**
 * \file StarState.h
 * \brief Provides the header for the StarState class.
 */

#ifndef StarState_cc
#define StarState_cc

#include "include/fe/GradientMatrix.h"
#include "include/riemann/StarState.h"

using namespace dealii;

/**
 * \brief Class for computing star state solutions in 1-D Riemann problems
 */
template <int dim>
class StarState
{
public:
  StarState(const std::shared_ptr<GradientMatrix<dim>> & gradient_matrix_,
            const DoFHandler<dim> & dof_handler_,
            const Triangulation<dim> & triangulation_,
            const unsigned int & n_components_);

  void compute_star_states(const Vector<double> & solution);

  void reinitialize();

  double get_star_state(const unsigned int & i, const unsigned int & j) const;

  virtual std::vector<double> compute(
    const std::vector<double> & solution_left,
    const std::vector<double> & solution_right,
    const Tensor<1, dim> & normal_vector) const = 0;

protected:
  /** \brief Pointer to gradient matrix object */
  std::shared_ptr<GradientMatrix<dim>> gradient_matrix;

  /** \brief Scalar 1st-order Lagrangian finite element */
  const FE_Q<dim> fe_scalar;

  /** \brief Pointer to cell quadrature */
  const DoFHandler<dim> * const dof_handler;

  /** \brief Degree of freedom handler for scalar case */
  DoFHandler<dim> dof_handler_scalar;

  /** \brief Number of solution components */
  unsigned int n_components;

  /** \brief Number of degrees of freedom for the vector case */
  unsigned int n_dofs;

  /** \brief Number of degrees of freedom for the scalar case */
  unsigned int n_dofs_scalar;

  /** \brief Sparsity pattern for vector system */
  SparsityPattern sparsity;

  /** \brief Sparsity pattern for scalar system */
  SparsityPattern sparsity_scalar;

  /** \brief Sparse matrix for star states */
  SparseMatrix<double> star_state_matrix;
};

#include "src/riemann/StarState.cc"

#endif
