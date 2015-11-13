/**
 * \file LaplacianDiffusion.h
 * \brief Provides the header for the LaplacianDiffusion class.
 */
#ifndef LaplacianDiffusion_h
#define LaplacianDiffusion_h

#include "Viscosity.h"

using namespace dealii;

/**
 * \class LaplacianDiffusion
 * \brief Class for applying Laplacian diffusion.
 */
template <int dim>
class LaplacianDiffusion : public ArtificialDiffusion<dim>
{
public:
  LaplacianDiffusion(
  const std::vector<FEValuesExtractors::Scalar> & scalar_extractors_,
  const std::vector<FEValuesExtractors::Vector> & vector_extractors_,
  const Viscosity<dim> * const viscosity,
  const unsigned int & n_quadrature_points_,
  const unsigned int & dofs_per_cell_);

  void apply(const Vector<double> & solution,
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

  /** \brief Viscosity */
  const Viscosity<dim> * const viscosity;

  /** \brief Number of quadrature points in cell */
  const unsigned int n_quadrature_points;

  /** \brief Number of degrees of freedom per cell */
  const unsigned int dofs_per_cell;
};

#include "LaplacianDiffusion.cc"

#endif
