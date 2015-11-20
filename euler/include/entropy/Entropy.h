/**
 * \file Entropy.h
 * \brief Provides the header for the Entropy class.
 */
#ifndef Entropy_h
#define Entropy_h

using namespace dealii;

#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

/**
 * \class Entropy
 * \brief Class for entropy.
 */
template <int dim>
class Entropy
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  Entropy(const bool & use_max_entropy_deviation_normalization,
          const double & domain_volume,
          const DoFHandler<dim> & dof_handler,
          const FESystem<dim> & fe,
          const QGauss<dim> & cell_quadrature,
          const QGauss<dim - 1> & face_quadrature);

  void reinitialize(const Vector<double> & solution);

  virtual std::vector<double> compute_entropy(
    const Vector<double> & solution,
    const FEValuesBase<dim> & fe_values) const = 0;

  virtual std::vector<double> compute_divergence_entropy_flux(
    const Cell & cell) = 0;

  virtual std::vector<Tensor<2, dim>> compute_entropy_flux_gradients_face(
    const Cell & cell, const unsigned int & i_face) = 0;

  virtual std::vector<Tensor<1, dim>> get_normal_vectors(
    const Cell & cell, const unsigned int & i_face) = 0;

  virtual std::vector<double> compute_entropy_normalization(
    const Vector<double> & solution, const Cell & cell) const = 0;

protected:
  std::vector<double> compute_max_entropy_deviation_normalization(
    const Vector<double> & solution, const Cell & cell) const;

  /** \brief Pointer to degree of freedom handler */
  const DoFHandler<dim> * const dof_handler;

  /** \brief Pointer to finite element system */
  const FESystem<dim> * const fe;

  /** \brief Pointer to cell quadrature */
  const QGauss<dim> * const cell_quadrature;

  /** \brief Pointer to face quadrature */
  const QGauss<dim - 1> * const face_quadrature;

  /** \brief Number of quadrature points per cell */
  const unsigned int n_q_points_cell;

  /** \brief Number of quadrature points per face */
  const unsigned int n_q_points_face;

private:
  virtual void reinitialize_group_fe_values(const Vector<double> & solution);

  void compute_average_entropy(const Vector<double> & solution);

  void compute_max_entropy_deviation(const Vector<double> & solution);

  /** \brief flag to use max entropy deviation from average as entropy
             normalization */
  const bool use_max_entropy_deviation_normalization;

  /** \brief Domain-average of entropy */
  double entropy_average;

  /** \brief Max entropy deviation from the average in the domain */
  double max_entropy_deviation;

  /** \brief Volume of the domain */
  const double domain_volume;
};

#include "src/entropy/Entropy.cc"

#endif
