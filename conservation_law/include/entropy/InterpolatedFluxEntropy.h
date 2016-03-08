/**
 * \file InterpolatedFluxEntropy.h
 * \brief Provides the header for the InterpolatedFluxEntropy class.
 */
#ifndef InterpolatedFluxEntropy_h
#define InterpolatedFluxEntropy_h

#include "include/entropy/Entropy.h"

using namespace dealii;

/**
 * \brief Abstract base class for entropy for which the entropy flux is
 *        interpolated.
 */
template <int dim>
class InterpolatedFluxEntropy : public Entropy<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Entropy<dim>::Cell;

  InterpolatedFluxEntropy(const bool & use_max_entropy_deviation_normalization_,
                          const double & domain_volume_,
                          const DoFHandler<dim> & dof_handler_,
                          const FESystem<dim> & fe_,
                          const QGauss<dim> & cell_quadrature_,
                          const QGauss<dim - 1> & face_quadrature_);

  std::vector<double> compute_entropy_residual(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const Cell & cell) override;

  double compute_max_entropy_jump(const Vector<double> & solution,
                                  const Cell & cell) override;

private:
  virtual std::vector<double> compute_entropy(
    const Vector<double> & solution,
    const FEValuesBase<dim> & fe_values) const = 0;

  virtual std::vector<double> compute_divergence_entropy_flux(
    const Cell & cell) = 0;

  virtual std::vector<Tensor<2, dim>> compute_entropy_flux_gradients_face(
    const Cell & cell, const unsigned int & i_face) = 0;

  virtual std::vector<Tensor<1, dim>> get_normal_vectors(
    const Cell & cell, const unsigned int & i_face) = 0;
};

#include "src/entropy/InterpolatedFluxEntropy.cc"

#endif
