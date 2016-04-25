/**
 * \file TransportAnalyticSolutionBounds.h
 * \brief Provides the header for the TransportAnalyticSolutionBounds class.
 */
#ifndef TransportAnalyticSolutionBounds_h
#define TransportAnalyticSolutionBounds_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>

#include "include/fct/FunctionDoFBounds.h"
#include "include/parameters/TransportProblemParameters.h"

using namespace dealii;

/**
 * \brief Abstract base class for implementing upper and lower bounds for
 *        the transport solution vector.
 */
template <int dim>
class TransportAnalyticSolutionBounds : public DoFBounds<dim>
{
public:
  TransportAnalyticSolutionBounds(
    TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const QGauss<dim> & cell_quadrature,
    const std::map<unsigned int, double> & dirichlet_values);

  void update(const Vector<double> & solution,
              const double & dt,
              const double & t_old);

protected:
  /** \brief cross section bounds */
  FunctionDoFBounds<dim> cross_section_bounds;

  /** \brief source bounds */
  FunctionDoFBounds<dim> source_bounds;

  /** \brief transport speed \f$v\f$ */
  const double speed;

  /** \brief map of DoF indices to Dirichlet BC values */
  const std::map<unsigned int, double> * const dirichlet_values;
};

#include "src/fct/TransportAnalyticSolutionBounds.cc"

#endif
