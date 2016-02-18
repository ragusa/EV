/**
 * \file ShallowWaterFCT.h
 * \brief Provides the header for the ShallowWaterFCT class.
 */

#ifndef ShallowWaterFCT_h
#define ShallowWaterFCT_h

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_accessor.h>
#include "include/parameters/RunParameters.h"
#include "include/postprocessing/PostProcessor.h"
#include "include/riemann/StarState.h"
#include "include/solvers/LinearSolver.h"

using namespace dealii;

/**
 * \brief Class for performing FCT for the shallow water equations.
 */
template <int dim>
class ShallowWaterFCT : public FCT<dim>
{
public:
  ShallowWaterFCT(const RunParameters<dim> & parameters,
                  const DoFHandler<dim> & dof_handler,
                  const Triangulation<dim> & triangulation,
                  const SparseMatrix<double> & lumped_mass_matrix,
                  const SparseMatrix<double> & consistent_mass_matrix,
                  const std::shared_ptr<StarState<dim>> & star_state,
                  const LinearSolver<dim> & linear_solver,
                  const SparsityPattern & sparsity_pattern,
                  const std::vector<unsigned int> & dirichlet_nodes,
                  const unsigned int & n_components,
                  const unsigned int & dofs_per_cell,
                  const std::vector<std::string> & component_names,
                  const bool & use_star_states_in_fct_bounds_,
                  const double & gravity);

protected:
  FullMatrix<double> compute_transformation_matrix(
    const Vector<double> & solution) const override;

  FullMatrix<double> compute_transformation_matrix_inverse(
    const Vector<double> & solution) const override;

  /** \brief Acceleration due to gravity */
  const double gravity;
};

#include "src/fct/ShallowWaterFCT.cc"
#endif
