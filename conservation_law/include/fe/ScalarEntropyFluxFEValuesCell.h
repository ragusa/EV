/**
 * \file ScalarEntropyFluxFEValuesCell.h
 * \brief Provides the header for the ScalarEntropyFluxFEValuesCell class.
 */
#ifndef ScalarEntropyFluxFEValuesCell_h
#define ScalarEntropyFluxFEValuesCell_h

#include <deal.II/base/quadrature_lib.h>
#include "include/fe/GroupFEValuesCell.h"

using namespace dealii;

/**
 * \class ScalarEntropyFluxFEValuesCell
 * \brief Class for computing entropy flux FE cell values for a scalar equation.
 */
template <int dim>
class ScalarEntropyFluxFEValuesCell : public GroupFEValuesCell<dim, false>
{
public:
  ScalarEntropyFluxFEValuesCell(const DoFHandler<dim> & solution_dof_handler,
                                const Triangulation<dim> & triangulation,
                                const QGauss<dim> & cell_quadrature)
    : GroupFEValuesCell<dim, false>(
        1, dim, solution_dof_handler, triangulation, cell_quadrature)
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
