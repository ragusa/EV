/**
 * \file ShallowWater2DDamBreakBC.h
 * \brief Provides the header for the ShallowWater2DDamBreakBC class.
 */
#ifndef ShallowWater2DDamBreakBC_h
#define ShallowWater2DDamBreakBC_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>
#include "include/bc/ShallowWaterBoundaryConditions.h"

using namespace dealii;

/**
 * \class ShallowWater2DDamBreakBC
 * \brief Class for 2-D dam break problem boundary conditions.
 */
template <int dim>
class ShallowWater2DDamBreakBC : public ShallowWaterBoundaryConditions<dim>
{
public:
  /** \brief Typedef for cell iterator */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  ShallowWater2DDamBreakBC(const FESystem<dim> & fe,
                           const QGauss<dim - 1> & face_quadrature,
                           const double & gravity);

private:
  void apply_boundary_condition(const Cell & cell,
                                const FEValues<dim> & fe_values_cell,
                                const FEFaceValues<dim> & fe_values_face,
                                const Vector<double> & solution,
                                const double & dt,
                                Vector<double> & cell_residual) override;

  void apply_interior_boundary_condition(const Cell & cell,
                                         const FEValues<dim> & fe_values_cell,
                                         const FEFaceValues<dim> & fe_values_face,
                                         const Vector<double> & solution,
                                         const double & dt,
                                         Vector<double> & cell_residual) override;
};

#include "src/bc/ShallowWater2DDamBreakBC.cc"

#endif
