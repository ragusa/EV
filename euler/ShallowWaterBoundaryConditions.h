/**
 * \file ShallowWaterBoundaryConditions.h
 * \brief Provides the header for the ShallowWaterBoundaryConditions class.
 */
#ifndef ShallowWaterBoundaryConditions_h
#define ShallowWaterBoundaryConditions_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

/**
 * \class ShallowWaterBoundaryConditions
 * \brief Base class for shallow water boundary conditions.
 */
template <int dim>
class ShallowWaterBoundaryConditions : BoundaryConditions<dim>
{
public:
  /** \brief Typedef for cell iterators */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  ShallowWaterBoundaryConditions(const FESystem<dim> & fe,
                                 const QGauss<dim> & face_quadrature,
                                 const double & gravity);

private:
  void integrate_face(const std::vector<double> & height,
                      const std::vector<Tensor<1, dim>> & momentum,
                      Vector<double> & cell_residual) const;

  /** \brief acceleration due to gravity */
  const double gravity;
};

#include "ShallowWaterBoundaryConditions.cc"

#endif
