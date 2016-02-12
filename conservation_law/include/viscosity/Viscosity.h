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
 *
 * The artificial viscosity \f$\nu\f$ derived from this class may be either the
 * type used by a standard Laplacian diffusion term:
 * \f[
 *   u_t + \nabla\cdot \mathbf{f}(u) - \nabla(\nu\nabla u) = 0 \,,
 * \f]
 * which makes the following contribution to the steady-state residual for
 * degree of freedom \f$i\f$, \f$r_i\f$:
 * \f[
 *   r_i = r_i + \int\limits_{S_i}\varphi_i\nabla(\nu\nabla u)dV \,,
 * \f]
 * or it may be the type used by a graph-theoretic diffusion term, which
 * makes the following contribution to the steady-state residual:
 * \f[
 *   r_i = r_i - \sum\limits_K\nu_K\sum\limits_{j\in J(K)}
 *     U_j b_K(\varphi_i, \varphi_j) \,,
 * \f]
 * where \f$b_K(\varphi_i, \varphi_j)\f$ is the local viscous bilinear form
 * for cell \f$K\f$. These viscosities differ by a factor of \f$h^2\f$; this
 * factor may be supplied as a viscosity multiplier if a viscosity is to
 * be used with its non-default diffusion type.
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
