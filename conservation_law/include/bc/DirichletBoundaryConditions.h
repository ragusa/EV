/**
 * \file DirichletBoundaryConditions.h
 * \brief Provides the header for the DirichletBoundaryConditions class.
 */
#ifndef DirichletBoundaryConditions_h
#define DirichletBoundaryConditions_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/vector.h>
#include "include/bc/BoundaryConditions.h"

using namespace dealii;

/**
 * \class DirichletBoundaryConditions
 * \brief Class for Dirichlet boundary conditions.
 *
 * This is a dummy class; nothing occurs when its apply() function is called.
 */
template <int dim>
class DirichletBoundaryConditions : public BoundaryConditions<dim>
{
public:
  /** \brief Typedef for cell iterator */
  using Cell = typename BoundaryConditions<dim>::Cell;

  DirichletBoundaryConditions(const FESystem<dim> & fe,
                              const QGauss<dim - 1> & face_quadrature)
    : BoundaryConditions<dim>(fe, face_quadrature)
  {
  }

private:
  void apply_boundary_condition(const Cell &,
                                const FEValues<dim> &,
                                const FEFaceValues<dim> &,
                                const Vector<double> &,
                                const double &,
                                Vector<double> &) override
  {
  }

  void apply_boundary_condition(const Cell &,
                                const FEValues<dim> &,
                                const FEFaceValues<dim> &,
                                const Vector<double> &,
                                const double &,
                                FullMatrix<double> &) override
  {
  }
};

#endif
