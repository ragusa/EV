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
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/lac/vector.h>

#include "include/base/ConservationLaw.h"
#include "include/entropy/TransportEntropy.h"
#include "include/fct/TransportExplicitEulerFCT.h"
#include "include/parameters/TransportRunParameters.h"
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
  Transport(const TransportRunParameters & params);

private:
  /** \brief Typedef for cell iterator */
  using Cell = typename ConservationLaw<dim>::Cell;

  /** \brief Typedef for face iterator */
  using Face = typename ConservationLaw<dim>::Face;

  std::vector<std::string> get_component_names() override;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_component_interpretations() override;

  void get_fe_extractors(
    std::vector<FEValuesExtractors::Scalar> & scalar_extractors,
    std::vector<FEValuesExtractors::Vector> & vector_extractors) const override;

  void define_problem() override;

  void assemble_lumped_mass_matrix() override;

  void perform_nonstandard_setup() override;

  void compute_inviscid_ss_matrix(const Vector<double> & solution,
                                  SparseMatrix<double> & matrix) override;

  void compute_ss_flux(const double & dt,
                       const Vector<double> & solution,
                       Vector<double> & ss_flux) override;

  void compute_ss_reaction(Vector<double> & ss_reaction) override;

  void compute_ss_rhs(const double & t, Vector<double> & ss_rhs) override;

  void update_flux_speeds() override;

  std::shared_ptr<Entropy<dim>> create_entropy() const override;

  std::shared_ptr<MaxWaveSpeed<dim>> create_max_wave_speed() const override;

  std::shared_ptr<ExplicitEulerFCT<dim>> create_fct() override;

  /** \brief Run parameters */
  TransportRunParameters transport_parameters;

  /** \brief Problem parameters */
  TransportProblemParameters<dim> problem_parameters;

  const FEValuesExtractors::Scalar extractor;
};

#include "src/base/Transport.cc"

#endif
