/**
 * \file ConstantViscosity.h
 * \brief Provides the header for the ConstantViscosity class.
 */
#ifndef ConstantViscosity_h
#define ConstantViscosity_h

#include "Viscosity.h"

using namespace dealii;

/**
 * \class ConstantViscosity
 * \brief Class for constant viscosity.
 */
template <int dim>
class ConstantViscosity : public Viscosity<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Viscosity<dim>::Cell;

  ConstantViscosity(const double & constant_value,
                    const DoFHandler<dim> & dof_handler);

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n) override;

private:
  const double constant_value;
};

#include "ConstantViscosity.cc"

#endif
