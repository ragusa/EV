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
template <int dim, typename FunctionType = double>
class GroupFEValuesCell : public GroupFEValuesBase<dim,FunctionType>
{
public:
  /**
   * \brief Typedef for cell iterators
   */
  using cell_iterator = typename GroupFEValuesBase<dim>::cell_iterator;

  GroupFEValuesCell(const unsigned int & n_components_solution,
                    const unsigned int & n_components_function,
                    const DoFHandler<dim> & solution_dof_handler,
                    const Triangulation<dim> & triangulation,
                    const QGauss<dim> & cell_quadrature,
                    const Vector<double> & solution,
                    const Vector<double> & aux_vector = Vector<double>());

  void reinit(const cell_iterator & solution_cell);

  void get_function_values(std::vector<double> & function_values) const;

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
