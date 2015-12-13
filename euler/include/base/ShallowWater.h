/**
 * \file ShallowWater.h
 * \brief Provides the header for the ShallowWater class.
 */
#ifndef ShallowWater_h
#define ShallowWater_h

#include <sys/stat.h>

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
#include "include/bc/ShallowWaterNoBC.h"
#include "include/bc/ShallowWaterSubcriticalOpenBC1D.h"
#include "include/bc/ShallowWaterSubcriticalWallBC1D.h"
#include "include/bc/ShallowWaterWallBC.h"
#include "include/entropy/Entropy.h"
#include "include/entropy/ShallowWaterEntropy.h"
#include "include/parameters/ShallowWaterParameters.h"
#include "include/parameters/ShallowWaterProblemParameters.h"
#include "include/postprocessing/ShallowWaterPostProcessor.h"
#include "include/postprocessing/ShallowWaterRiemannSolver.h"
#include "include/viscosity/MaxWaveSpeed.h"
#include "include/viscosity/ShallowWaterMaxWaveSpeed.h"
#include "include/viscosity/ViscosityMultiplier.h"
#include "include/viscosity/ShallowWaterViscosityMultiplier.h"

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

  void assemble_lumped_mass_matrix() override;

  void define_problem() override;

  void perform_nonstandard_setup() override;

  void compute_ss_residual(const double & dt, Vector<double> & solution) override;

  void update_flux_speeds() override;

  void compute_inviscid_fluxes(const std::vector<double> & height,
                               const std::vector<Tensor<1, dim>> & momentum,
                               std::vector<Tensor<1, dim>> & height_flux,
                               std::vector<Tensor<2, dim>> & momentum_flux) const;

  std::vector<Tensor<1, dim>> compute_velocity(
    const std::vector<double> & height,
    const std::vector<Tensor<1, dim>> & momentum) const;

  std::vector<double> compute_speed(
    const std::vector<double> & height,
    const std::vector<Tensor<1, dim>> & momentum) const;

  std::vector<double> compute_sound_speed(
    const std::vector<double> & height) const;

  std::shared_ptr<Entropy<dim>> create_entropy() const override;

  std::shared_ptr<MaxWaveSpeed<dim>> create_max_wave_speed() const override;

  std::shared_ptr<ViscosityMultiplier<dim>> create_viscosity_multiplier()
    const override;

  std::shared_ptr<DataPostprocessor<dim>> create_auxiliary_postprocessor()
    const override;

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

#include "src/base/ShallowWater.cc"

#endif
