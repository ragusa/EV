/**
 * \file HighOrderViscosity.h
 * \brief Provides the header for the HighOrderViscosity class.
 */
#ifndef HighOrderViscosity_h
#define HighOrderViscosity_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include "Viscosity.h"

using namespace dealii;

/**
 * \class HighOrderViscosity
 * \brief Class for high-order viscosity
 */
template <int dim>
class HighOrderViscosity : public Viscosity<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Viscosity<dim>::Cell;

  HighOrderViscosity(const std::shared_ptr<Viscosity<dim>> & low_order_viscosity,
                     const std::shared_ptr<Viscosity<dim>> & entropy_viscosity,
                     const bool & use_low_order_for_first_step,
                     const DoFHandler<dim> & dof_handler);

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n) override;

private:
  /** \brief Pointer to low-order viscosity */
  std::shared_ptr<Viscosity<dim>> low_order_viscosity;

  /** \brief Pointer to entropy viscosity */
  std::shared_ptr<Viscosity<dim>> entropy_viscosity;

  /** \brief Flag to use low-order viscosity for first time step */
  const bool use_low_order_for_first_step;
};

#include "HighOrderViscosity.cc"

#endif
