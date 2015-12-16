/**
 * \file BoundaryConditions.h
 * \brief Provides the header for the BoundaryConditions class.
 */
#ifndef BoundaryConditions_h
#define BoundaryConditions_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>

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

  /** \brief Typedef for face iterators */
  typedef typename DoFHandler<dim>::face_iterator Face;

  BoundaryConditions(const FESystem<dim> & fe,
                     const QGauss<dim - 1> & face_quadrature);

  void apply(const Cell & cell,
             const FEValues<dim> & fe_values,
             const Vector<double> & solution,
             const double & dt,
             Vector<double> & cell_residual);

  virtual void apply_dirichlet_boundary_conditions() {}

protected:
  virtual void apply_boundary_condition(const Cell & cell,
                                        const FEValues<dim> & fe_values_cell,
                                        const FEFaceValues<dim> & fe_values_face,
                                        const Vector<double> & solution,
                                        const double & dt,
                                        Vector<double> & cell_residual) = 0;

  /** \brief Finite element system */
  const FESystem<dim> fe;

  /** \brief Quadrature for face */
  const QGauss<dim - 1> face_quadrature;

  /** \brief FE values for face */
  FEFaceValues<dim> fe_values_face;

  /** \brief Number of faces per cell */
  const unsigned int faces_per_cell;

  /** \brief Number of degrees of freedom per cell */
  const unsigned int dofs_per_cell;

  /** \brief Number of quadrature points for face */
  const unsigned int n_quadrature_points_face;
};

#include "src/bc/BoundaryConditions.cc"

#endif
