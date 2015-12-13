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
 * This class adds a standard Laplacian diffusion term:
 * \f[
 *   u_t + \nabla\cdot \mathbf{f}(u) - \nabla(\nu\nabla u) = 0 \,,
 * \f]
 * which makes the following contribution to the steady-state residual for
 * degree of freedom \f$i\f$, \f$r_i\f$:
 * \f[
 *   r_i = r_i + \int\limits_{S_i}\varphi_i\nabla(\nu\nabla u)dV \,,
 * \f]
 */
template <int dim>
class LaplacianDiffusion : public ArtificialDiffusion<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename ArtificialDiffusion<dim>::Cell;

  LaplacianDiffusion(
    const std::vector<FEValuesExtractors::Scalar> & scalar_extractors_,
    const std::vector<FEValuesExtractors::Vector> & vector_extractors_,
    const unsigned int & n_quadrature_points_,
    const unsigned int & dofs_per_cell_);

  void apply(std::shared_ptr<Viscosity<dim>> viscosity,
             const Vector<double> & solution,
             const Cell & cell,
             const FEValues<dim> & fe_values,
             Vector<double> & cell_residual) const override;

private:
  /** \brief Vector of scalar extractors for conservation law */
  const std::vector<FEValuesExtractors::Scalar> scalar_extractors;

  /** \brief Vector of vector extractors for conservation law */
  const std::vector<FEValuesExtractors::Vector> vector_extractors;

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
