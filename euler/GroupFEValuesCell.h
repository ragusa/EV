/**
 * \file GroupFEValuesCell.h
 * \brief Provides the header for the GroupFEValuesCell class.
 */
#ifndef GroupFEValuesCell_h
#define GroupFEValuesCell_h

#include <deal.II/base/quadrature_lib.h>
#include "GroupFEValuesBase.h"

using namespace dealii;

/**
 * \class GroupFEValuesCell
 * \brief Class for computing group FEM cell values for a function.
 */
template <int dim, bool is_scalar = true>
class GroupFEValuesCell : public GroupFEValuesBase<dim, is_scalar>
{
public:
  /** \brief Typedef for cell iterators */
  using Cell = typename GroupFEValuesBase<dim, is_scalar>::Cell;

  /** \brief Typedef for function value type (scalar or vector) */
  using ValueType = typename GroupFEValuesBase<dim, is_scalar>::ValueType;

  /** \brief Typedef for function gradient type (1st or 2nd-order tensor) */
  using GradientType = typename GroupFEValuesBase<dim, is_scalar>::GradientType;

  GroupFEValuesCell(const unsigned int & n_components_solution,
                    const unsigned int & n_components_function,
                    const DoFHandler<dim> & solution_dof_handler,
                    const Triangulation<dim> & triangulation,
                    const QGauss<dim> & cell_quadrature,
                    const Vector<double> & aux_vector = Vector<double>());

  void reinit(const Cell & solution_cell);

  void get_function_values(std::vector<ValueType> & function_values) const;

  void get_function_gradients(
    std::vector<GradientType> & function_gradients) const;

  std::vector<double> get_function_divergences() const;

protected:
  /** \brief Cell quadrature */
  const QGauss<dim> cell_quadrature;

  /** \brief Number of quadrature points */
  const unsigned int n_quadrature_points;

  /** \brief Finite Element cell values for function */
  FEValues<dim> function_fe_values_cell;
};

#include "GroupFEValuesCell.cc"

#endif
