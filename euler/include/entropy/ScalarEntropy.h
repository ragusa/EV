/**
 * \file ScalarEntropy.h
 * \brief Provides the header for the ScalarEntropy class.
 */
#ifndef ScalarEntropy_h
#define ScalarEntropy_h

#include "include/entropy/Entropy.h"
#include "include/fe/ScalarEntropyFluxFEValuesCell.h"
#include "include/fe/ScalarEntropyFluxFEValuesFace.h"

using namespace dealii;

/**
 * \class ScalarEntropy
 * \brief Class for entropy for a scalar conservation law equation.
 */
template <int dim>
class ScalarEntropy : public Entropy<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Entropy<dim>::Cell;

  ScalarEntropy(const double & domain_volume,
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
  ScalarEntropyFluxFEValuesCell<dim> entropy_flux_fe_values_cell;

  ScalarEntropyFluxFEValuesFace<dim> entropy_flux_fe_values_face;
};

#include "src/entropy/ScalarEntropy.cc"

#endif
