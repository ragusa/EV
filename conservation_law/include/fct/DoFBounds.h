/**
 * \file DoFBounds.h
 * \brief Provides the header for the DoFBounds class.
 */
#ifndef DoFBounds_h
#define DoFBounds_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/vector.h>

using namespace dealii;

/**
 * \brief Abstract base class for implementing upper and lower bounds for
 *        some degree of freedom indexed vector such as a solution vector.
 */
template <int dim>
class DoFBounds
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  DoFBounds(const DoFHandler<dim> & dof_handler, const FESystem<dim> & fe);

  void compute_min_max_dof_vector(const Vector<double> & dof_vector,
                                  Vector<double> & min_dof_vector,
                                  Vector<double> & max_dof_vector) const;

  void widen(const DoFBounds<dim> & other_dof_bounds);

  void shrink(const DoFBounds<dim> & other_dof_bounds);

  void remove_bounds(const std::vector<double> & dof_indices);

  /** \brief lower bound vector */
  Vector<double> lower;

  /** \brief upper bound vector */
  Vector<double> upper;

protected:
  /** \brief degree of freedom handler */
  const DoFHandler<dim> * const dof_handler;

  /** \brief finite element system */
  const FESystem<dim> * const fe;

  /** \brief number of degrees of freedom */
  const unsigned int n_dofs;

  /** \brief number of degrees of freedom per cell */
  const unsigned int dofs_per_cell;
};

#include "src/fct/DoFBounds.cc"

#endif
