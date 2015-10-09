/**
 * \file ShallowWaterEntropyFluxFEValuesFace.h
 * \brief Provides the header for the ShallowWaterEntropyFluxFEValuesFace class.
 */
#ifndef ShallowWaterEntropyFluxFEValuesFace_h
#define ShallowWaterEntropyFluxFEValuesFace_h

#include <deal.II/base/quadrature_lib.h>
#include "GroupFEValuesFace.h"

using namespace dealii;

/**
 * \class ShallowWaterEntropyFluxFEValuesFace
 * \brief Class for computing entropy flux FE face values for the shallow
 *        water equations.
 */
template <int dim>
class ShallowWaterEntropyFluxFEValuesFace : public GroupFEValuesFace<dim>
{
public:
  ShallowWaterEntropyFluxFEValuesFace(const unsigned int & n_components,
                    const DoFHandler<dim> & solution_dof_handler,
                    const Triangulation<dim> & triangulation,
                    const QGauss<dim> & face_quadrature,
                    const Vector<double> & solution,
                    const Vector<double> & bathymetry_vector)
  {
    this->compute_function_dof_values();
  }

protected:
  double function(const std::vector<double> & solution, const double & bathymetry) const override
  {
    double entropy_flux = 1.0;
    for (unsigned int i = 0; i < solution.size(); ++i)
      value *= solution[i];
    return value;
  }
};

#include "ShallowWaterEntropyFluxFEValuesFace.cc"

#endif
