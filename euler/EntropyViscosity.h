/**
 * \file EntropyViscosity.h
 * \brief Provides the header for the EntropyViscosity class.
 */
#ifndef EntropyViscosity_h
#define EntropyViscosity_h

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include "Entropy.h"
#include "Viscosity.h"
#include "GroupFEValuesCell.h"

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

  EntropyViscosity(const DoFHandler<dim> & dof_handler);

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n) override;

private:
  std::shared_ptr<Entropy<dim>> entropy;

  const double residual_coefficient;

  const double jump_coefficient;

  const CellMap * const cell_diameter;

  const unsigned int n_q_points_cell;

  const FESystem<dim> * fe;

  const QGauss<dim> * cell_quadrature;
};

#include "EntropyViscosity.cc"

#endif
