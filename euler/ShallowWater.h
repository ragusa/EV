/**
 * \file ShallowWater.h
 * \brief Provides the header for the ShallowWater class.
 */
#ifndef ShallowWater_h
#define ShallowWater_h

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
#include "ShallowWaterEntropyFluxFEValuesFace.h"
#include "ShallowWaterParameters.h"
#include "ShallowWaterPostProcessor.h"
#include "ShallowWaterRiemannSolver.h"

/**
 * \brief Class for solving the shallow water equations.
 *
 * The shallow water equations are the following:
 * \f[
 *   \frac{\partial h}{\partial t}
 *   + \nabla\cdot\mathbf{q} = 0 ,
 * \f]
 * \f[
 *   \frac{\partial\mathbf{q}}{\partial t}
 *   + \nabla\cdot\left(\mathbf{q}\otimes\mathbf{v}
 *   + \frac{1}{2}g h^2 \mathbf{I}\right) =
 *   - g h \nabla b ,
 * \f]
 * where \f$h\f$ is water height, \f$\mathbf{v}\f$ is velocity,
 * \f$\mathbf{q}=h\mathbf{v}\f$ is the ``momentum'',
 * \f$g\f$ is acceleration due to gravity, \f$\mathbf{I}\f$ is the identity
 * tensor, and \f$b\f$ is the bathymetry (bottom topography).
 */
template <int dim>
class ShallowWater : public ConservationLaw<dim>
{
public:
  ShallowWater(const ShallowWaterParameters<dim> & params);

private:
  /**
   * \brief Typedef for cell iterators
   */
  using cell_iterator = typename ConservationLaw<dim>::cell_iterator;

  std::vector<std::string> get_component_names() override;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_component_interpretations() override;

  void assemble_lumped_mass_matrix() override;

  void define_problem() override;

  void perform_additional_setup() override;

  void compute_ss_residual(Vector<double> & solution) override;

  void update_flux_speeds() override;

  void compute_entropy(const Vector<double> & solution,
                       const FEValuesBase<dim> & fe_values,
                       Vector<double> & entropy) const override;

  void compute_divergence_entropy_flux(
    const Vector<double> & solution,
    const FEValuesBase<dim> & fe_values,
    Vector<double> & divergence_entropy_flux) const override;

  void update_entropy_viscosities(const double & dt) override;

  double compute_max_entropy_jump(const Vector<double> & solution,
                                          const cell_iterator & cell) const;

  void compute_inviscid_fluxes(const std::vector<double> & height,
                               const std::vector<Tensor<1, dim>> & momentum,
                               std::vector<Tensor<1, dim>> & height_flux,
                               std::vector<Tensor<2, dim>> & momentum_flux) const;

  void compute_viscous_fluxes(
    const double & viscosity,
    const std::vector<Tensor<1, dim>> & height_gradient,
    const std::vector<Tensor<2, dim>> & momentum_gradient,
    std::vector<Tensor<1, dim>> & height_viscous_flux,
    std::vector<Tensor<2, dim>> & momentum_viscous_flux) const;

  void compute_velocity(const std::vector<double> & height,
                        const std::vector<Tensor<1, dim>> & momentum,
                        std::vector<Tensor<1, dim>> & velocity) const;

  void compute_entropy_derivative_height(
    const std::vector<double> & height,
    const std::vector<Tensor<1, dim>> & momentum,
    std::vector<double> & entropy_derivative_height) const;

  void compute_entropy_derivative_momentum(
    const std::vector<double> & height,
    const std::vector<Tensor<1, dim>> & momentum,
    std::vector<Tensor<1,dim>> & entropy_derivative_momentum) const;

  void output_results(PostProcessor<dim> & postprocessor) const override;

  /** \brief Parameters for shallow water equations */
  ShallowWaterParameters<dim> sw_parameters;

  /** \brief FE values extractor for height */
  const FEValuesExtractors::Scalar height_extractor;

  /** \brief FE values extractor for momentum */
  const FEValuesExtractors::Vector momentum_extractor;

  /** \brief Acceleration due to gravity \f$g\f$ */
  double gravity;

  /** \brief Finite element for bathymetry */
  const FE_Q<dim> fe_bathymetry;

  /** \brief Degree of freedom handler for bathymetry */
  DoFHandler<dim> dof_handler_bathymetry;

  /** \brief Bathymetry (bottom topography) vector \f$b\f$ */
  Vector<double> bathymetry_vector;

  /** \brief Bathymetry (bottom topography) function \f$b\f$ */
  std::shared_ptr<Function<dim>> bathymetry_function;
};

#include "ShallowWater.cc"

#endif
