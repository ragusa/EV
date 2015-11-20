/**
 * \file EntropyViscosity.h
 * \brief Provides the header for the EntropyViscosity class.
 */
#ifndef EntropyViscosity_h
#define EntropyViscosity_h

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include "include/parameters/ConservationLawParameters.h"
#include "include/entropy/Entropy.h"
#include "include/viscosity/Viscosity.h"
#include "include/fe/GroupFEValuesCell.h"

using namespace dealii;

/**
 * \class EntropyViscosity
 * \brief Class for entropy viscosity.
 */
template <int dim>
class EntropyViscosity : public Viscosity<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Viscosity<dim>::Cell;

  /** \brief Alias for cell iterator map to double */
  using CellMap = typename Viscosity<dim>::CellMap;

  EntropyViscosity(const ConservationLawParameters<dim> & parameters,
                   const std::shared_ptr<Entropy<dim>> & entropy,
                   const CellMap & cell_diameter,
                   const FESystem<dim> & fe,
                   const DoFHandler<dim> & dof_handler,
                   const QGauss<dim> & cell_quadrature,
                   const QGauss<dim - 1> & face_quadrature);

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n) override;

private:
  std::vector<double> compute_entropy_residual(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const Cell & cell) const;

  double compute_max_entropy_jump(const Cell & cell) const;

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

  /** \brief Pointer to map of cell iterator to cell diameter */
  const CellMap * const cell_diameter;

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
};

#include "src/viscosity/EntropyViscosity.cc"

#endif
