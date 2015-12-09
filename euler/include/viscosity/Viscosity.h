/**
 * \file Viscosity.h
 * \brief Provides the header for the Viscosity class.
 */
#ifndef Viscosity_h
#define Viscosity_h

#include <map>

#include <deal.II/dofs/dof_handler.h>

using namespace dealii;

/**
 * \class Viscosity
 * \brief Class for artificial viscosity.
 */
template <int dim>
class Viscosity
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  /** \brief Alias for cell iterator map to double */
  using CellMap = std::map<Cell, double>;

  Viscosity(const DoFHandler<dim> & dof_handler);

  virtual void update(const Vector<double> & new_solution,
                      const Vector<double> & old_solution,
                      const double & dt,
                      const unsigned int & n) = 0;

  virtual void reinitialize();

  double & operator[](const Cell & cell);

  unsigned int size() const;

  CellMap get_values() const;

protected:
  /** \brief Map of cell iterator to viscosity value */
  CellMap values;

  /** \brief Pointer to degree of freedom handler */
  const DoFHandler<dim> * const dof_handler;
};

#include "src/viscosity/Viscosity.cc"

#endif