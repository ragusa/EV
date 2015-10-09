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
class ShallowWaterEntropyFluxFEValuesFace : public GroupFEValuesFace<dim, false>
{
public:
  ShallowWaterEntropyFluxFEValuesFace(
    const DoFHandler<dim> & solution_dof_handler,
    const Triangulation<dim> & triangulation,
    const QGauss<dim-1> & face_quadrature,
    const Vector<double> & solution,
    const Vector<double> & bathymetry_vector,
    const double & gravity)
    : GroupFEValuesFace<dim, false>(dim + 1,
                                    dim,
                                    solution_dof_handler,
                                    triangulation,
                                    face_quadrature,
                                    solution,
                                    bathymetry_vector),
      gravity(gravity)
  {
    this->compute_function_dof_values();
  }

protected:
  /**
   * \brief Computes entropy flux at a single point.
   *
   * \param[in] solution solution at one point
   * \param[in] bathymetry bathymetry function \f$b\f$ at one point
   *
   * \return entropy flux vector
   */
  std::vector<double> function(const std::vector<double> & solution,
                               const double & bathymetry) const override
  {
    // extract solution components from solution vector
    const double height = solution[0];
    Tensor<1, dim> momentum;
    for (unsigned int d = 0; d < dim; ++d)
      momentum[d] = solution[d + 1];

    // compute dot product of momentum
    const double momentum_dot_product = momentum * momentum;

    // compute each component (dimension) of entropy flux
    std::vector<double> entropy_flux(dim);
    for (unsigned int d = 0; d < dim; ++d)
      entropy_flux[d] = gravity * (height + bathymetry) * momentum[d] +
        0.5 * momentum_dot_product / height / height * momentum[d];

    return entropy_flux;
  }

  /** acceleration due to gravity \f$g\f$ */
  const double gravity;
};

#endif
