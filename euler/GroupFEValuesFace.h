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
template <int dim, bool is_scalar = true>
class GroupFEValuesFace : public GroupFEValuesBase<dim, is_scalar>
{
public:
  /** \brief Typedef for cell iterators */
  using Cell = typename GroupFEValuesBase<dim, is_scalar>::Cell;

  /** \brief Typedef for function value type (scalar or vector) */
  using ValueType = typename GroupFEValuesBase<dim, is_scalar>::ValueType;

  /** \brief Typedef for function gradient type (1st or 2nd-order tensor) */
  using GradientType = typename GroupFEValuesBase<dim, is_scalar>::GradientType;

  GroupFEValuesFace(const unsigned int & n_components_solution,
                    const unsigned int & n_components_function,
                    const DoFHandler<dim> & solution_dof_handler,
                    const Triangulation<dim> & triangulation,
                    const QGauss<dim - 1> & face_quadrature,
                    const Vector<double> & aux_vector = Vector<double>());

  void reinit(const Cell & solution_cell, const unsigned int & face);

  void get_function_values(std::vector<ValueType> & function_values) const;

  void get_function_gradients(
    std::vector<GradientType> & function_gradients) const;

  std::vector<Tensor<1, dim>> get_normal_vectors() const;

protected:
  /** \brief Face quadrature */
  const QGauss<dim - 1> face_quadrature;

  /** \brief Number of quadrature points */
  const unsigned int n_quadrature_points;

  /** \brief Finite Element face values for function */
  FEFaceValues<dim> function_fe_values_face;
};

#include "GroupFEValuesFace.cc"

#endif
