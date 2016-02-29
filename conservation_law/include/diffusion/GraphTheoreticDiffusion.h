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
 * This class computes a graph-theoretic diffusion matrix:
 * \f[
 *   D_{i,j} = \sum\limits_{K\subset S_{i,j}} \nu_K b_K(\varphi_i, \varphi_j) \,,
 * \f]
 * where \f$\nu_K\f$ is the cell viscosity, and
 * \f$b_K(\varphi_i,\varphi_j)\f$ is the local viscous bilinear form:
 * \f[
 *   b_K(\varphi_i,\varphi_j) = \left\{\begin{array}{c l}
 *     |K|                & j = i, \quad i,j\in\mathcal{I}(K) \\
 *     -\frac{|K|}{n_K-1} & j \ne i, \quad i,j\in\mathcal{I}(K), \quad m(j)=m(i)
 * \\
 *     0 & \mbox{otherwise}
 *   \end{array}\right. \,,
 * \f]
 * where \f$|K|\f$ is the cell volume, \f$n_K=\mbox{card}(\mathcal{I}(K))\f$,
 * and \f$m(i)\f$ is the component index of degree of freedom \f$i\f$.
 */
template <int dim>
class GraphTheoreticDiffusion : public ArtificialDiffusion<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  GraphTheoreticDiffusion(const DoFHandler<dim> & dof_handler_,
                          const unsigned int & dofs_per_cell,
                          const unsigned int & n_components);

  void compute_diffusion_matrix(const Vector<double> & solution,
                                const std::shared_ptr<Viscosity<dim>> viscosity,
                                SparseMatrix<double> & diffusion_matrix) override;

protected:
  /** \brief Pointer to degree of freedom handler */
  const DoFHandler<dim> * const dof_handler;

  /** \brief Number of degrees of freedom per cell */
  const unsigned int dofs_per_cell;

  /** \brief Number of solution components */
  const unsigned int n_components;

  /** \brief Number of degrees of freedom per cell per solution component */
  const unsigned int dofs_per_cell_per_component;

  /** \brief Unit cell matrix for local bilinear forms, \f$\mathrm{B}^K\f$*/
  FullMatrix<double> unit_cell_matrix;
};

#include "src/diffusion/GraphTheoreticDiffusion.cc"

#endif
