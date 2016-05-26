/**
 * \file TransportUpwindAnalyticSolutionBounds.h
 * \brief Provides the header for the TransportUpwindAnalyticSolutionBounds class.
 */
#ifndef TransportUpwindAnalyticSolutionBounds_h
#define TransportUpwindAnalyticSolutionBounds_h

#include <algorithm>
#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>

#include "include/geometry/BoundaryDistance.h"
#include "include/parameters/TransportProblemParameters.h"

using namespace dealii;

/**
 * \brief Class for implementing upwind transport solution bounds.
 */
template <int dim>
class TransportUpwindAnalyticSolutionBounds : public DoFBounds<dim>
{
public:
  TransportUpwindAnalyticSolutionBounds(
    TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const unsigned int & n_samples);

  void update(const Vector<double> & solution,
              const double & s,
              const double & t_old);

protected:
  /** \brief cross section function */
  const FunctionParser<dim> * const cross_section_function;

  /** \brief source function */
  FunctionParser<dim> * const source_function;

  /** \brief incoming solution value (assumed constant) */
  const double incoming;

  /** \brief transport direction */
  const Tensor<1, dim> direction;

  /** \brief positions of each DoF */
  std::vector<Point<dim>> x;

  /** \brief boundary distance function */
  const std::shared_ptr<BoundaryDistance<dim>> boundary_distance;

  /** \brief number of sampling points for min/max determination */
  const unsigned int n_samples;

  /** \brief sample point vector */
  std::vector<Point<dim>> x_sample;

  /** \brief fractions for sampling */
  std::vector<double> alpha_sample;

  /** \brief cross section samples */
  std::vector<double> sigma_samples;

  /** \brief source samples */
  std::vector<double> source_samples;
};

#include "src/fct/TransportUpwindAnalyticSolutionBounds.cc"

#endif
