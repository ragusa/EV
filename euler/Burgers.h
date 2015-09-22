/**
 * \file Burgers.h
 * \brief Provides the header for the Burgers class.
 */
#ifndef Burgers_h
#define Burgers_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/lac/vector.h>
#include "ConservationLaw.h"
#include "BurgersParameters.h"

/**
 * \brief Class for solving the Burgers equation.
 *
 * The Burgers equation is the following:
 * \f[
 *   \frac{\partial u}{\partial t}
 *   + \nabla\cdot\left(\frac{1}{2}u^2\mathbf{v}\right) = 0 ,
 * \f]
 * or
 * \f[
 *   \frac{\partial u}{\partial t}
 *   + u\mathbf{v}\cdot\nabla u = 0 ,
 * \f]
 * where \f$\mathbf{v}\f$ is a constant velocity field:
 * \f[
 *   \mathbf{v} = \left[\begin{array}{c}1\end{array}\right]
 *   \qquad \mbox{(1-D)}
 * \f]
 * \f[
 *   \mathbf{v} = \left[\begin{array}{c}1\\1\end{array}\right]
 *   \qquad \mbox{(2-D)}
 * \f]
 * \f[
 *   \mathbf{v} = \left[\begin{array}{c}1\\1\\1\end{array}\right]
 *   \qquad \mbox{(3-D)}
 * \f]
 */
template <int dim>
class Burgers : public ConservationLaw<dim>
{
public:

    Burgers(const BurgersParameters<dim> &params);

private:

    BurgersParameters<dim> burgers_parameters;

    const FEValuesExtractors::Scalar velocity_extractor;

    std::vector<std::string> get_component_names() override;

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
       get_component_interpretations() override;

    void assemble_lumped_mass_matrix() override;

    void define_problem() override;

    void compute_ss_residual(Vector<double> &solution) override;

/*
    void compute_face_ss_residual(
      FEFaceValues<dim> &fe_face_values,
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      Vector<double> &cell_residual) override;
*/

    //void compute_ss_jacobian() override;

    void update_flux_speeds() override;

    void compute_entropy(
      const Vector<double>    &solution,
      const FEValuesBase<dim> &fe_values,
      Vector<double>          &entropy) const override;

    void compute_divergence_entropy_flux(
      const Vector<double> &solution,
      const FEValues<dim>  &fe_values,
      Vector<double>       &divergence_entropy_flux) const override;
};

#include "Burgers.cc"

#endif
