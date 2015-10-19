/**
 * \file ShallowWaterSubcriticalOpenBC1D.h
 * \brief Provides the header for the ShallowWaterSubcriticalOpenBC1D class.
 */
#ifndef ShallowWaterSubcriticalOpenBC1D_h
#define ShallowWaterSubcriticalOpenBC1D_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>
#include "ShallowWaterBoundaryConditions.h"

using namespace dealii;

/**
 * \class ShallowWaterSubcriticalOpenBC1D
 * \brief Class for subcritical open boundary conditions in 1-D for shallow water
 *        equations
 */
template <int dim>
class ShallowWaterSubcriticalOpenBC1D : public BoundaryConditions<dim>
{
public:
  /** \brief Typedef for cell iterator */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  ShallowWaterSubcriticalOpenBC1D(const FESystem<dim> & fe,
                                  const QGauss<dim - 1> & face_quadrature,
                                  const double & gravity,
                                  const double & height_left,
                                  const double & height_right);

private:
  void apply_boundary_condition(const Cell & cell,
                                const FEValues<dim> & fe_values_cell,
                                const FEFaceValues<dim> & fe_values_face,
                                const Vector<double> & solution,
                                Vector<double> & cell_residual) override;

  /** \brief left boundary value for height */
  const double height_left;

  /** \brief right boundary value for height */
  const double height_right;
};

#include "ShallowWaterSubcriticalOpenBC1D.cc"

#endif
