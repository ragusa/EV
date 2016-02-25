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
#include "include/viscosity/ConstantMaxWaveSpeed.h"

/**
 * \brief Class for solving a scalar transport equation.
 *
 * The scalar transport equation is the following:
 * \f[
 *   \frac{1}{v}\frac{\partial u}{\partial t}
 *   + \nabla\cdot\left(\mathbf{u\Omega}\right)
 *   + \sigma u = q \,,
 * \f]
 * where \f$u\f$ is the conserved quantity, \f$v\f$ is the speed of the
 * waves or particles, \f$\mathbf{\Omega}\f$ is the transport direction,
 * \f$\sigma\f$ is the reaction coefficient (cross section), and \f$q\f$
 * is a positive source term.
 */
template <int dim>
class Transport : public ConservationLaw<dim>
{
public:
  Transport(const TransportParameters<dim> & params);

private:
  std::vector<std::string> get_component_names() override;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_component_interpretations() override;

  void get_fe_extractors(
    std::vector<FEValuesExtractors::Scalar> & scalar_extractors,
    std::vector<FEValuesExtractors::Vector> & vector_extractors) const override;

  void assemble_lumped_mass_matrix() override;

  void perform_nonstandard_setup() override;

  //void set_boundary_ids();

  void define_problem() override;

  void compute_ss_flux(const double & dt,
                       const Vector<double> & solution,
                       Vector<double> & ss_flux) override;

  void compute_ss_reaction(Vector<double> & ss_reaction) override;

  void compute_ss_rhs(const double & t, Vector<double> & ss_rhs) override;

  void update_flux_speeds() override;

  std::shared_ptr<Entropy<dim>> create_entropy() const override;

  std::shared_ptr<MaxWaveSpeed<dim>> create_max_wave_speed() const override;

  TransportParameters<dim> transport_parameters;

  const FEValuesExtractors::Scalar extractor;

  /** \brief Transport speed \f$v\f$ */
  //double transport_speed;

  /** \brief x-component of transport direction, \f$\Omega_x\f$ */
  //double transport_direction_x;

  /** \brief y-component of transport direction, \f$\Omega_y\f$ */
  //double transport_direction_y;

  /** \brief z-component of transport direction, \f$\Omega_z\f$ */
  //double transport_direction_z;

  /** \brief Transport direction \f$\mathbf{\Omega}\f$ */
  //Tensor<1, dim> transport_direction;

  /** \brief Function parser of cross section \f$\sigma\f$ */
  //FunctionParser<dim> cross_section_function;

  /** \brief Function parser of source \f$q\f$ */
  //FunctionParser<dim> source_function;
};

#include "src/base/Transport.cc"

#endif
