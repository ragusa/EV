/**
 * \file GroupFEValues.h
 * \brief Provides the header for the GroupFEValues class.
 */
#ifndef GroupFEValues_h
#define GroupFEValues_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

/**
 * \class GroupFEValues
 * \brief Class for computing group FEM values for a function.
 *
 * It is assumed that both the solution and the function use linear Langrangian
 * finite elements, i.e., FE_Q<dim> elements of degree 1.
 */
template <int dim>
class GroupFEValues
{
public:
  GroupFEValues(
    const unsigned int & n_components,
    const Vector<double> & solution,
    const Triangulation<dim> & triangulation,
    const QGauss<dim> & cell_quadrature);

  std::vector<double> get_function_values(
    const typename DoFHandler<dim>::active_cell_iterator & cell);

  std::vector<double> get_function_dof_values() const;

protected:
  /**
   * \brief Evaluates the function at a DoF support point.
   *
   * This function must be overridden by derived classes.
   *
   * \param[in] solution all components of solution at a DoF support point
   *
   * \return function value at a DoF support point
   */
  virtual double function(const std::vector<double> & solution) const = 0;

  void compute_function_dof_values();

  /** \brief Exception for a modulus not returning zero */
  DeclException2(ExcModulusNotZero, int, int,
    << arg1 << " % " << arg2 << " != 0");

  /** \brief Number of solution components */
  const unsigned int n_components;

  /** \brief Finite element for function */
  const FE_Q<dim> fe;

  /** \brief DoF handler for function */
  DoFHandler<dim> dof_handler;

  /** \brief Number of function DoFs per cell */
  const unsigned int n_function_dofs_per_cell;

  /** \brief Number of solution DoFs per cell */
  const unsigned int n_solution_dofs_per_cell;

  /** \brief Cell quadrature */
  const QGauss<dim> cell_quadrature;

  /** \brief Number of quadrature points */
  const unsigned int n_quadrature_points;

  /** \brief Cell iterator for function DoF handler */
  typename DoFHandler<dim>::active_cell_iterator cell_function;

  /** \brief End cell iterator for function DoF handler */
  typename DoFHandler<dim>::active_cell_iterator endc_function;

  /** \brief Number of DoFs for function */
  unsigned int n_dofs;

  /** \brief Function values at DoF support points */
  std::vector<double> function_dof_values;

  /** \brief Pointer to solution vector */
  const Vector<double> * const solution;

  /** \brief Finite Element values for function */
  FEValues<dim> fe_values_function;
};

#include "GroupFEValues.cc"

#endif
