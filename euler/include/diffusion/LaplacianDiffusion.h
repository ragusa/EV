/**
 * \file LaplacianDiffusion.h
 * \brief Provides the header for the LaplacianDiffusion class.
 */
#ifndef LaplacianDiffusion_h
#define LaplacianDiffusion_h

#include "include/viscosity/Viscosity.h"

using namespace dealii;

/**
 * \class LaplacianDiffusion
 * \brief Class for applying Laplacian diffusion.
 *
 * This class computes a diffusion matrix representing a standard Laplacian
 * diffusion term:
 * \f[
 *   u_t + \nabla\cdot \mathbf{f}(u) - \nabla(\nu\nabla u) = 0 \,.
 * \f]
 * This diffusion matrix is computed as
 * \f[
 *   D_{i,j} = -\sum\limits_{K\subset S_{i,j}}
 *     \int\limits_K\varphi_i(\mathbf{x})\nabla\cdot(\nu_K\nabla\varphi_j(\mathbf{x}))
 * dV \,,
 * \f]
 * which after integration by parts and dropping the boundary term, is
 * \f[
 *   D_{i,j} = \sum\limits_{K\subset S_{i,j}}
 *     \nu_K\int\limits_K\nabla\varphi_i(\mathbf{x})\cdot\nabla\varphi_j(\mathbf{x})
 * dV \,,
 * \f]
 * where \f$\nu_K\f$ is the cell viscosity.
 */
template <int dim>
class LaplacianDiffusion : public ArtificialDiffusion<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename ArtificialDiffusion<dim>::Cell;

  LaplacianDiffusion(
    const std::vector<FEValuesExtractors::Scalar> & scalar_extractors,
    const std::vector<FEValuesExtractors::Vector> & vector_extractors,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const QGauss<dim> & cell_quadrature,
    const unsigned int & dofs_per_cell);

  void compute_diffusion_matrix(const Vector<double> & solution,
                                const std::shared_ptr<Viscosity<dim>> viscosity,
                                SparseMatrix<double> & diffusion_matrix) override;

private:
  /** \brief Vector of scalar extractors for conservation law */
  const std::vector<FEValuesExtractors::Scalar> scalar_extractors;

  /** \brief Vector of vector extractors for conservation law */
  const std::vector<FEValuesExtractors::Vector> vector_extractors;

  /** \brief Pointer to degree of freedom handler */
  const DoFHandler<dim> * const dof_handler;

  /** \brief Pointer to finite element system */
  const FESystem<dim> * const fe;

  /** \brief Pointer to cell quadrature */
  const QGauss<dim> * const cell_quadrature;

  /** \brief Number of scalar components in solution */
  const unsigned int n_scalar_components;

  /** \brief Number of vector components in solution */
  const unsigned int n_vector_components;

  /** \brief Number of quadrature points in cell */
  const unsigned int n_quadrature_points;

  /** \brief Number of degrees of freedom per cell */
  const unsigned int dofs_per_cell;
};

#include "src/diffusion/LaplacianDiffusion.cc"

#endif
