/**
 * \file Transport.h
 * \brief Provides the header for the Transport class.
 */
#ifndef Transport_h
#define Transport_h

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
#include "include/parameters/TransportParameters.h"
#include "include/parameters/TransportProblemParameters.h"
#include "include/viscosity/TransportMaxWaveSpeed.h"

/**
 * \brief Class for solving a scalar transport equation.
 *
 * The scalar transport equation is the following:
 * \f[
 *   \frac{\partial u}{\partial t}
 *   + \nabla\cdot\left(\mathbf{v}u\right)
 *   + \sigma u = q ,
 * \f]
 */
template <int dim>
class Transport : public ConservationLaw<dim>
{
public:
  Transport(const TransportParameters<dim> & params);

private:
  TransportParameters<dim> transport_parameters;

  const FEValuesExtractors::Scalar extractor;

  std::vector<std::string> get_component_names() override;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_component_interpretations() override;

  void get_fe_extractors(
    std::vector<FEValuesExtractors::Scalar> & scalar_extractors,
    std::vector<FEValuesExtractors::Vector> & vector_extractors) const override;

  void assemble_lumped_mass_matrix() override;

  void define_problem() override;

  void compute_ss_residual(const double & dt, Vector<double> & solution) override;

  void update_flux_speeds() override;

  std::shared_ptr<Entropy<dim>> create_entropy() const override;

  std::shared_ptr<MaxWaveSpeed<dim>> create_max_wave_speed() const override;
};

#include "src/base/Transport.cc"

#endif
