/**
 * \file ScalarEntropyFluxFEValuesFace.h
 * \brief Provides the header for the ScalarEntropyFluxFEValuesFace class.
 */
#ifndef ScalarEntropyFluxFEValuesFace_h
#define ScalarEntropyFluxFEValuesFace_h

#include <deal.II/base/quadrature_lib.h>
#include "include/fe/GroupFEValuesFace.h"

using namespace dealii;

/**
 * \class ScalarEntropyFluxFEValuesFace
 * \brief Class for computing entropy flux FE face values for a scalar equation.
 */
template <int dim>
class ScalarEntropyFluxFEValuesFace : public GroupFEValuesFace<dim, false>
{
public:
  ScalarEntropyFluxFEValuesFace(const DoFHandler<dim> & solution_dof_handler,
                                const Triangulation<dim> & triangulation,
                                const QGauss<dim - 1> & face_quadrature)
    : GroupFEValuesFace<dim, false>(
        1, dim, solution_dof_handler, triangulation, face_quadrature)
  {
  }

protected:
  /**
   * \brief Computes entropy flux at a single point.
   *
   * \param[in] solution solution at one point
   * \param[in] aux auxiliary parameter, not used here
   *
   * \return entropy flux vector
   */
  std::vector<double> function(const std::vector<double> & solution,
                               const double &) const override
  {
    // compute each component (dimension) of entropy flux
    std::vector<double> entropy_flux(dim);
    for (unsigned int d = 0; d < dim; ++d)
      entropy_flux[d] = solution[0];

    return entropy_flux;
  }
};

#endif
