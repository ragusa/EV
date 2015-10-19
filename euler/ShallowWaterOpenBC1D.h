/**
 * \file ShallowWaterOpenBC1D.h
 * \brief Provides the header for the ShallowWaterOpenBC1D class.
 */
#ifndef ShallowWaterOpenBC1D_h
#define ShallowWaterOpenBC1D_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

/**
 * \class ShallowWaterOpenBC1D
 * \brief Class for open boundary conditions in 1-D for shallow water
 *        equations
 */
template <int dim>
class ShallowWaterOpenBC1D : BoundaryConditions<dim>
{
public:
  /** \brief Typedef for cell iterators */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  ShallowWaterOpenBC1D(FESystem<dim> & fe,
                       const QGauss<dim> & face_quadrature,
                       const double & gravity);

private:
  void apply_boundary_condition(const Cell & cell,
                                const FEValues<dim> & fe_values_cell,
                                const FEValuesFace<dim> & fe_values_face,
                                const Vector<double> & solution,
                                Vector<double> & cell_residual) override;
};

#include "ShallowWaterOpenBC1D.cc"

#endif
