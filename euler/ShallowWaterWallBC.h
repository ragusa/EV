/**
 * \file ShallowWaterWallBC.h
 * \brief Provides the header for the ShallowWaterWallBC class.
 */
#ifndef ShallowWaterWallBC_h
#define ShallowWaterWallBC_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>
#include "ShallowWaterBoundaryConditions.h"

using namespace dealii;

/**
 * \class ShallowWaterWallBC
 * \brief Class for wall boundary conditions for shallow water equations.
 */
template <int dim>
class ShallowWaterWallBC : public ShallowWaterBoundaryConditions<dim>
{
public:
  /** \brief Typedef for cell iterator */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  ShallowWaterWallBC(const FESystem<dim> & fe,
                     const QGauss<dim - 1> & face_quadrature,
                     const double & gravity);

private:
  void apply_boundary_condition(const Cell & cell,
                                const FEValues<dim> & fe_values_cell,
                                const FEFaceValues<dim> & fe_values_face,
                                const Vector<double> & solution,
                                const double & dt,
                                Vector<double> & cell_residual) override;
};

#include "ShallowWaterWallBC.cc"

#endif
