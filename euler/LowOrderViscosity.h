/**
 * \file LowOrderViscosity.h
 * \brief Provides the header for the LowOrderViscosity class.
 */
#ifndef LowOrderViscosity_h
#define LowOrderViscosity_h

using namespace dealii;

/**
 * \class LowOrderViscosity
 * \brief Class for low-order viscosity.
 *
 * The low-order viscosity is computed as
 * \f[
 *   \nu_K^L = c_{max} h_K \lambda_{K,max} \,,
 * \f]
 * where \f$\lambda_{K,max}\f$ is the maximum flux speed on cell \f$K\f$.
 */
template <int dim>
class LowOrderViscosity : public Viscosity<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Viscosity<dim>::Cell;

  /** \brief Alias for cell iterator map to double */
  using CellMap = typename Viscosity<dim>::CellMap;

  LowOrderViscosity(const double & c_max,
                    CellMap & cell_diameter,
                    CellMap & max_flux_speed,
                    const DoFHandler<dim> & dof_handler);

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n) override;

private:
  const double c_max;

  CellMap * const cell_diameter;

  CellMap * const max_flux_speed;
};

#include "LowOrderViscosity.cc"

#endif
