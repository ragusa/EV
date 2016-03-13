/**
 * \file EntropyViscosity.h
 * \brief Provides the header for the EntropyViscosity class.
 */
#ifndef EntropyViscosity_h
#define EntropyViscosity_h

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include "include/parameters/RunParameters.h"
#include "include/entropy/Entropy.h"
#include "include/viscosity/Viscosity.h"
#include "include/fe/GroupFEValuesCell.h"

using namespace dealii;

/**
 * \brief Class for entropy viscosity.
 *
 * When using Laplacian diffusion, the entropy viscosity for cell \f$K\f$
 * at time \f$n\f$ is computed as
 * \f[
 *   \nu^{\eta,n}_K \equiv h_K^2 \max\limits_{q\in Q(K)}\left(
 *     \frac{c_R |R_q^n| + c_J J_K^n}{\hat{\eta}_q}
 *     \right)\,.
 * \f]
 * When using graph-theoretic diffusion, the leading \f$h_K^2\f$ coefficient
 * is dropped:
 * \f[
 *   \nu^{\eta,n}_K \equiv \max\limits_{q\in Q(K)}\left(
 *     \frac{c_R |R_q^n| + c_J J_K^n}{\hat{\eta}_q}
 *     \right)\,.
 * \f]
 */
template <int dim>
class EntropyViscosity : public Viscosity<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Viscosity<dim>::Cell;

  /** \brief Alias for cell iterator map to double */
  using CellMap = typename Viscosity<dim>::CellMap;

  // EntropyViscosity(const RunParameters & parameters,
  EntropyViscosity(const RunParameters & parameters,
                   const std::shared_ptr<Entropy<dim>> & entropy,
                   const FESystem<dim> & fe,
                   const DoFHandler<dim> & dof_handler,
                   const QGauss<dim> & cell_quadrature,
                   const QGauss<dim - 1> & face_quadrature,
                   const bool & use_in_laplacian_term);

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n) override;

private:
  double compute_viscosity_multiplier_laplacian(const Cell & cell) const;

  double compute_viscosity_multiplier_graphtheoretic(const Cell & cell) const;

  void smooth_entropy_viscosity_max();

  void smooth_entropy_viscosity_average();

  /** \brief Pointer to entropy */
  std::shared_ptr<Entropy<dim>> entropy;

  /** \brief Coefficient for entropy residual */
  const double residual_coefficient;

  /** \brief Coefficient for entropy flux jumps */
  const double jump_coefficient;

  /** \brief Weighting to be used if an average weighting is to be applied */
  const double smoothing_weight;

  /** \brief Pointer to finite element system */
  const FESystem<dim> * fe;

  /** \brief Pointer to cell quadrature */
  const QGauss<dim> * cell_quadrature;

  /** \brief Pointer to face quadrature */
  const QGauss<dim - 1> * face_quadrature;

  /** \brief Number of quadrature points per cell */
  const unsigned int n_q_points_cell;

  /** \brief Number of quadrature points per face */
  const unsigned int n_q_points_face;

  /** \brief Number of faces per cell */
  const unsigned int faces_per_cell;

  /** \brief Function pointer for computing viscosity multiplier */
  double (EntropyViscosity<dim>::*compute_viscosity_multiplier_ptr)(
    const Cell & cell) const;
};

#include "src/viscosity/EntropyViscosity.cc"

#endif
