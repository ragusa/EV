/**
 * \file BoundaryConditions.h
 * \brief Provides the header for the BoundaryConditions class.
 */
#ifndef BoundaryConditions_h
#define BoundaryConditions_h

#include <deal.II/dofs/dof_handler.h>
//#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
/*
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/grid/tria.h>
*/
#include <deal.II/lac/vector.h>
//#include "Exceptions.h"

using namespace dealii;

/**
 * \class BoundaryConditions
 * \brief Base class for boundary conditions.
 */
template <int dim>
class BoundaryConditions
{
public:
  /** \brief Typedef for cell iterators */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  BoundaryConditions();

  void apply(const Cell & cell, const FEValues<dim> & fe_values,
    const Vector<double> & solution);

private:
  virtual void apply_boundary_condition(const Cell & cell,
    const FEValues<dim> & fe_values_cell, const FEValuesFace<dim> & fe_values_face,
    const Vector<double> & solution, Vector<double> & cell_residual) = 0;
  
  /** \brief Finite element system for function */
  const FESystem<dim> fe;
};

#include "BoundaryConditions.cc"

#endif
