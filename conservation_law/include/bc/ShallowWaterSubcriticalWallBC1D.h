/**
 * \file ShallowWaterSubcriticalWallBC1D.h
 * \brief Provides the header for the ShallowWaterSubcriticalWallBC1D class.
 */
#ifndef ShallowWaterSubcriticalWallBC1D_h
#define ShallowWaterSubcriticalWallBC1D_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>
#include "include/bc/ShallowWaterBoundaryConditions.h"

using namespace dealii;

/**
 * \class ShallowWaterSubcriticalWallBC1D
 * \brief Class for subcritical wall boundary conditions in 1-D for shallow water
 *        equations
 */
template <int dim>
class ShallowWaterSubcriticalWallBC1D : public ShallowWaterBoundaryConditions<dim>
{
public:
  /** \brief Typedef for cell iterator */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  ShallowWaterSubcriticalWallBC1D(const FESystem<dim> & fe,
                                  const QGauss<dim - 1> & face_quadrature,
                                  const double & gravity);

private:
  void apply_boundary_condition(const Cell & cell,
                                const FEValues<dim> & fe_values_cell,
                                const FEFaceValues<dim> & fe_values_face,
                                const Vector<double> & solution,
                                const double & dt,
                                Vector<double> & cell_residual) override;

  void get_interior_values(const double & face_position,
                           const double & wave_velocity,
                           const double & dt,
                           const Cell & cell,
                           const Vector<double> & solution,
                           double & speed_of_sound_interior,
                           double & velocity_interior) const;
};

#include "src/bc/ShallowWaterSubcriticalWallBC1D.cc"

#endif
