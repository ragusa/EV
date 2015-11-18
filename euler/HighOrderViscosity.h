/**
 * \file HighOrderViscosity.h
 * \brief Provides the header for the HighOrderViscosity class.
 */
#ifndef HighOrderViscosity_h
#define HighOrderViscosity_h

using namespace dealii;

/**
 * \class HighOrderViscosity
 * \brief Class for high-order viscosity
 */
template <int dim>
class HighOrderViscosity
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Viscosity<dim>::Cell;

  HighOrderViscosity(const std::shared_ptr<Viscosity<dim>> & low_order_viscosity,
                     const std::shared_ptr<Viscosity<dim>> & entropy_viscosity,
                     const bool & use_low_order_for_first_step,
                     const DoFHandler<dim> & dof_handler);

private:
  std::shared_ptr<Viscosity<dim>> low_order_viscosity;

  std::shared_ptr<Viscosity<dim>> entropy_viscosity;

  const bool use_low_order_for_first_step;
};

#include "HighOrderViscosity.cc"

#endif
