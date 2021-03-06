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
#include "include/base/ConservationLaw.h"
#include "include/entropy/ScalarEntropy.h"
#include "include/parameters/BurgersRunParameters.h"
#include "include/parameters/BurgersProblemParameters.h"
#include "include/viscosity/BurgersMaxWaveSpeed.h"

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
  Burgers(const BurgersRunParameters & params);

private:
  BurgersRunParameters burgers_parameters;

  const FEValuesExtractors::Scalar velocity_extractor;

  std::vector<std::string> get_component_names() override;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_component_interpretations() override;

  void get_fe_extractors(
    std::vector<FEValuesExtractors::Scalar> & scalar_extractors,
    std::vector<FEValuesExtractors::Vector> & vector_extractors) const override;

  void assemble_lumped_mass_matrix() override;

  void define_problem() override;

  void compute_ss_flux(const double & dt,
                       const Vector<double> & solution,
                       Vector<double> & ss_flux) override;

  void compute_ss_rhs(const double & t, Vector<double> & ss_rhs) override;

  void update_flux_speeds() override;

  std::shared_ptr<Entropy<dim>> create_entropy() const override;

  std::shared_ptr<MaxWaveSpeed<dim>> create_max_wave_speed() const override;
};

#include "src/base/Burgers.cc"

#endif
