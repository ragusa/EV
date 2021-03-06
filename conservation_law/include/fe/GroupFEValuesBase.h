/**
 * \file GroupFEValuesBase.h
 * \brief Provides the header for the GroupFEValuesBase class.
 */
#ifndef GroupFEValuesBase_h
#define GroupFEValuesBase_h

#include <deal.II/base/exceptions.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>
#include "include/other/Exceptions.h"

using namespace dealii;

/**
 * \class GroupFEValuesBase
 * \brief Base class for computing group FEM values for a function.
 *
 * It is assumed that both the solution and the function use linear Langrangian
 * finite elements, i.e., FE_Q<dim> elements of degree 1.
 */
template <int dim, bool is_scalar = true>
class GroupFEValuesBase
{
public:
  /** \brief Typedef for cell iterators */
  typedef typename DoFHandler<dim>::active_cell_iterator Cell;

  /**
   * \brief Typedef for function FE values extractor.
   *
   * Depending on whether the function returns a scalar, i.e., a double,
   * or a dim-dimensional vector, i.e., a Tensor<1,dim>, the FE values
   * extractor must be of type FEValuesExtractors::Scalar or
   * FEValuesExtractors::Vector, respectively.
   */
  typedef
    typename std::conditional<is_scalar,
                              FEValuesExtractors::Scalar,
                              FEValuesExtractors::Vector>::type ExtractorType;

  /** \brief Typedef for function value type (scalar or vector) */
  typedef
    typename std::conditional<is_scalar, double, Tensor<1, dim>>::type ValueType;

  /** \brief Typedef for function gradient type (1st or 2nd-order tensor) */
  typedef
    typename std::conditional<is_scalar, Tensor<1, dim>, Tensor<2, dim>>::type
      GradientType;

  GroupFEValuesBase(const unsigned int & n_components_solution,
                    const unsigned int & n_components_function,
                    const DoFHandler<dim> & solution_dof_handler,
                    const Triangulation<dim> & triangulation,
                    const Vector<double> & aux_vector = Vector<double>());

  Vector<double> get_function_dof_values() const;

  void reinitialize(const Vector<double> & solution);

protected:
  /**
   * \brief Evaluates the function at a DoF support point.
   *
   * This function must be overridden by derived classes.
   *
   * \param[in] solution all components of solution at a DoF support point
   * \param[in] aux optional auxiliary quantity at a DoF support point
   *
   * \return function value at a DoF support point
   */
  virtual std::vector<double> function(const std::vector<double> & solution,
                                       const double & aux = 0.0) const = 0;

  /** \brief Number of solution components */
  const unsigned int n_components_solution;

  /** \brief Number of function components */
  const unsigned int n_components_function;

  /** \brief Finite element system for function */
  const FESystem<dim> fe;

  /** \brief FE values extractor for function */
  const ExtractorType function_extractor;

  /** \brief DoF handler for function */
  DoFHandler<dim> dof_handler;

  /** \brief Number of function DoFs per cell */
  const unsigned int n_function_dofs_per_cell;

  /** \brief Number of solution DoFs per cell */
  const unsigned int n_solution_dofs_per_cell;

  /** \brief Map from solution cell iterator to function cell iterator */
  std::map<Cell, Cell> solution_cell_to_function_cell_map;

  /** \brief Number of DoFs for function */
  unsigned int n_dofs;

  /** \brief Number of support points */
  unsigned int n_support_points;

  /** \brief Function values at DoF support points */
  Vector<double> function_dof_values;

  /** \brief Pointer to optional auxiliary vector */
  const Vector<double> * const aux_vector;

  /** \brief Flag to signal that function DoF values have been initialized */
  bool function_values_initialized;
};

#include "src/fe/GroupFEValuesBase.cc"

#endif
