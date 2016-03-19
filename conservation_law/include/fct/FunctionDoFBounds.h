/**
 * \file FunctionDoFBounds.h
 * \brief Provides the header for the FunctionDoFBounds class.
 */
#ifndef FunctionDoFBounds_h
#define FunctionDoFBounds_h

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include "include/fct/DoFBounds.h"

using namespace dealii;

/**
 * \brief Abstract base class for implementing upper and lower bounds for
 *        a function at degree of freedom indexed points.
 */
template <int dim>
class FunctionDoFBounds : public DoFBounds<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFBounds<dim>::Cell;

  FunctionDoFBounds(Function<dim> & function,
                    const bool & function_is_time_dependent,
                    const DoFHandler<dim> & dof_handler,
                    const FESystem<dim> & fe,
                    const QGauss<dim> & cell_quadrature);

  void update(const double & time);

protected:
  void compute_function_bounds();

  /** \brief Function \f$f(\mathbf{x},t)\f$ */
  Function<dim> * const function;

  /** \brief Flag that function is time-dependent */
  const bool function_is_time_dependent;

  /** \brief cell quadrature */
  const QGauss<dim> * const cell_quadrature;

  /** \brief number of quadrature points per cell */
  const unsigned int n_q_points_cell;
};

#include "src/fct/FunctionDoFBounds.cc"

#endif
