/**
 * \file GraphTheoreticDiffusion.h
 * \brief Provides the header for the GraphTheoreticDiffusion class.
 */
#ifndef GraphTheoreticDiffusion_h
#define GraphTheoreticDiffusion_h

#include "include/diffusion/ArtificialDiffusion.h"
#include "include/viscosity/Viscosity.h"

using namespace dealii;

/**
 * \class GraphTheoreticDiffusion
 * \brief Class for graph-theoretic diffusion.
 *
 * This class adds a graph-theoretic diffusion term, which
 * makes the following contribution to the steady-state residual:
 * \f[
 *   r_i = r_i - \sum\limits_K\nu_K\sum\limits_{j\in J(K)}
 *     U_j b_K(\varphi_i, \varphi_j) \,,
 * \f]
 * where \f$b_K(\varphi_i, \varphi_j)\f$ is the local viscous bilinear form
 * for cell \f$K\f$.
 */
template <int dim>
class GraphTheoreticDiffusion : public ArtificialDiffusion<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  GraphTheoreticDiffusion(const unsigned int & dofs_per_cell,
                          const unsigned int & n_components);

  void apply(std::shared_ptr<Viscosity<dim>> viscosity,
             const Vector<double> & solution,
             const Cell & cell,
             const FEValues<dim> & fe_values,
             Vector<double> & cell_residual) const override;

protected:
  /** \brief Number of degrees of freedom per cell */
  const unsigned int dofs_per_cell;

  /** \brief Number of solution components */
  const unsigned int n_components;

  /** \brief Number of degrees of freedom per cell per solution component */
  const unsigned int dofs_per_cell_per_component;

  /** \brief Sparsity pattern for the viscous sums matrix */
  SparsityPattern sparsity;

  /** \brief Matrix of sums of graph-theoretic viscous bilinear forms */
  SparseMatrix<double> viscous_sums;
};

#include "src/diffusion/GraphTheoreticDiffusion.cc"

#endif
