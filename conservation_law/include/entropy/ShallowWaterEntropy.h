/**
 * \file ShallowWaterEntropy.h
 * \brief Provides the header for the ShallowWaterEntropy class.
 */
#ifndef ShallowWaterEntropy_h
#define ShallowWaterEntropy_h

#include <string>
#include "include/parameters/ShallowWaterRunParameters.h"
#include "include/fe/ShallowWaterEntropyFluxFEValuesCell.h"
#include "include/fe/ShallowWaterEntropyFluxFEValuesFace.h"

using namespace dealii;

/**
 * \class ShallowWaterEntropy
 * \brief Class for entropy for the shallow water equations.
 */
template <int dim>
class ShallowWaterEntropy : public Entropy<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Entropy<dim>::Cell;

  ShallowWaterEntropy(const ShallowWaterRunParameters<dim> & parameters,
                      const FEValuesExtractors::Scalar & height_extractor,
                      const FEValuesExtractors::Vector & momentum_extractor,
                      const double & gravity,
                      const Vector<double> & bathymetry_vector,
                      const double & domain_volume,
                      const DoFHandler<dim> & dof_handler,
                      const FESystem<dim> & fe,
                      const Triangulation<dim> & triangulation,
                      const QGauss<dim> & cell_quadrature,
                      const QGauss<dim - 1> & face_quadrature);

  std::vector<double> compute_entropy(
    const Vector<double> & solution,
    const FEValuesBase<dim> & fe_values) const override;

  std::vector<double> compute_divergence_entropy_flux(const Cell & cell) override;

  std::vector<Tensor<2, dim>> compute_entropy_flux_gradients_face(
    const Cell & cell, const unsigned int & i_face) override;

  std::vector<Tensor<1, dim>> get_normal_vectors(
    const Cell & cell, const unsigned int & i_face) override;

  std::vector<double> compute_entropy_normalization(
    const Vector<double> & solution, const Cell & cell) const override;

  void reinitialize_group_fe_values(const Vector<double> & solution) override;

private:
  std::vector<double> compute_local_entropy_normalization(
    const Vector<double> & solution, const Cell & cell) const;

  /** \brief Function pointer for computing entropy normalization */
  std::vector<double> (
    ShallowWaterEntropy<dim>::*compute_entropy_normalization_ptr)(
    const Vector<double> & solution, const Cell & cell) const;

  const FEValuesExtractors::Scalar height_extractor;

  const FEValuesExtractors::Vector momentum_extractor;

  const double gravity;

  ShallowWaterEntropyFluxFEValuesCell<dim> entropy_flux_fe_values_cell;

  ShallowWaterEntropyFluxFEValuesFace<dim> entropy_flux_fe_values_face;

  const std::string entropy_normalization_option;
};

#include "src/entropy/ShallowWaterEntropy.cc"

#endif
