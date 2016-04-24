/**
 * \file CharacteristicFCTFilter.h
 * \brief Provides the header for the CharacteristicFCTFilter class.
 */
#ifndef CharacteristicFCTFilter_h
#define CharacteristicFCTFilter_h

#include "include/fct/ExplicitEulerFCTFilter.h"

using namespace dealii;

/**
 * \brief Abstract base class for characteristic FCT filters.
 */
template <int dim>
class CharacteristicFCTFilter : public ExplicitEulerFCTFilter<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  CharacteristicFCTFilter(
    const RunParameters & run_parameters,
    const std::shared_ptr<Limiter<dim>> limiter,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const SparseMatrix<double> & lumped_mass_matrix,
    const std::map<unsigned int, double> & dirichlet_values_);

  virtual void filter_antidiffusive_fluxes(
    const Vector<double> & old_solution,
    const double & dt,
    const Vector<double> & inviscid_ss_flux,
    const Vector<double> & ss_reaction,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const Vector<double> & ss_rhs,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix) override;

  virtual bool check_bounds(const Vector<double> & new_solution) override;

protected:
  virtual void compute_solution_bounds(const Vector<double> & old_solution,
                                       const double & dt,
                                       const Vector<double> & ss_reaction,
                                       const Vector<double> & ss_rhs) override;

  virtual void compute_antidiffusion_bounds(
    const Vector<double> & old_solution,
    const double & dt,
    const Vector<double> & inviscid_ss_flux,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const Vector<double> & ss_rhs) override;

  void transform_vector(const Vector<double> & solution,
                        const Vector<double> & vector_original,
                        Vector<double> & vector_transformed) const;

  void transform_matrix(const Vector<double> & solution,
                        const SparseMatrix<double> & matrix_original,
                        SparseMatrix<double> & matrix_transformed) const;

  void create_dof_indices_lists();

  virtual FullMatrix<double> compute_transformation_matrix(
    const Vector<double> & solution) const = 0;

  /** \brief sparsity pattern */
  SparsityPattern sparsity_pattern;

  /** \brief matrix of transformed antidiffusive fluxes \f$\hat{\mathbf{P}}\f$ */
  SparseMatrix<double> characteristic_antidiffusion_matrix;

  /** \brief old solution vector in characteristic variables \hat{\mathbf{U}}^n */
  Vector<double> old_solution_characteristic;

  /** \brief temporary vector */
  Vector<double> tmp_vector;

  /** \brief list of node index of each degree of freedom */
  std::vector<unsigned int> node_indices;

  /** \brief list of component index of each degree of freedom */
  std::vector<unsigned int> component_indices;

  /** \brief list of degree of freedom indices for each node */
  std::vector<std::vector<unsigned int>> node_dof_indices;

  /** \brief number of nodes */
  unsigned int n_nodes;
};

#include "src/fct/CharacteristicFCTFilter.cc"

#endif
