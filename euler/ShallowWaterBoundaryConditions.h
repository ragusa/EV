/**
 * \file ShallowWaterBoundaryConditions.h
 * \brief Provides the header for the ShallowWaterBoundaryConditions class.
 */
#ifndef ShallowWaterBoundaryConditions_h
#define ShallowWaterBoundaryConditions_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

/**
 * \class ShallowWaterBoundaryConditions
 * \brief Base class for shallow water boundary conditions.
 */
template <int dim>
class ShallowWaterBoundaryConditions
{
public:
  /** \brief Typedef for cell iterators */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  ShallowWaterBoundaryConditions();

  void apply(const Cell & cell, const FEValues<dim> & fe_values,
    const Vector<double> & solution);

private:
  virtual void apply_boundary_condition(const Cell & cell,
    const FEValues<dim> & fe_values_cell, const FEValuesFace<dim> & fe_values_face,
    const Vector<double> & solution, Vector<double> & cell_residual) const = 0;

  void integrate_face() const;
};

#include "ShallowWaterBoundaryConditions.cc"

#endif
