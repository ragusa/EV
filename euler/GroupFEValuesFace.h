/**
 * \file GroupFEValuesFace.h
 * \brief Provides the header for the GroupFEValuesFace class.
 */
#ifndef GroupFEValuesFace_h
#define GroupFEValuesFace_h

#include <deal.II/base/quadrature_lib.h>
#include "GroupFEValuesBase.h"

using namespace dealii;

/**
 * \class GroupFEValuesFace
 * \brief Class for computing group FEM face values for a function.
 */
template <int dim>
class GroupFEValuesFace : public GroupFEValuesBase<dim>
{
public:
  /**
   * \brief Typedef for cell iterators
   */
  using cell_iterator = typename GroupFEValuesBase<dim>::cell_iterator;

  GroupFEValuesFace(const unsigned int & n_components,
                    const DoFHandler<dim> & solution_dof_handler,
                    const Triangulation<dim> & triangulation,
                    const QGauss<dim> & face_quadrature,
                    const Vector<double> & solution,
                    const Vector<double> & aux_vector = Vector<double>());

  void reinit(const cell_iterator & solution_cell, const unsigned int & face);

  void get_function_values(std::vector<double> & function_values) const;

protected:
  /** \brief Face quadrature */
  const QGauss<dim> face_quadrature;

  /** \brief Number of quadrature points */
  const unsigned int n_quadrature_points;

  /** \brief Finite Element face values for function */
  FEFaceValues<dim> function_fe_values_face;
};

#include "GroupFEValuesFace.cc"

#endif
