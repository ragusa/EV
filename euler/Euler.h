/** \file Euler.h
 *  \brief Provides the header for the Euler class.
 */
#ifndef Euler_h
#define Euler_h

#include <deal.II/base/symmetric_tensor.h>

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
#include "EulerParameters.h"

/** \class Euler
 *  \brief Provides the necessary functions for Euler equations.
 */
template <int dim>
class Euler : public ConservationLaw<dim>
{
public:

  Euler(const EulerParameters<dim> &params);

private:

  EulerParameters<dim> euler_parameters;

  // number of components and position of components in solution vector
  static const unsigned int n_euler_components = dim+2;

  const FEValuesExtractors::Scalar density_extractor;
  const FEValuesExtractors::Vector momentum_extractor;
  const FEValuesExtractors::Scalar energy_extractor;

  // vector of names of each component
  std::vector<std::string> get_component_names();

  // data component interpretation (scalar or vector) for outputting solution
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
     get_component_interpretations();

  void define_problem();
  void output_solution(double time);

  void compute_ss_residual(Vector<double> &solution) override;
  void compute_cell_ss_residual(FEValues<dim> &fe_values,
                                FEFaceValues<dim> &fe_face_values,
                                const typename DoFHandler<dim>::active_cell_iterator &cell,
                                Vector<double> &cell_residual);
/*
  void compute_face_ss_residual(FEFaceValues<dim> &fe_face_values,
                                const typename DoFHandler<dim>::active_cell_iterator &cell,
                                Vector<double> &cell_residual);
*/
  void compute_ss_jacobian();
  
  void compute_inviscid_fluxes(
    const std::vector<double>         & density,
    const std::vector<Tensor<1,dim> > & momentum,
    const std::vector<double>         & energy,
    std::vector<Tensor<1,dim> > & density_flux,
    std::vector<Tensor<2,dim> > & momentum_flux,
    std::vector<Tensor<1,dim> > & energy_flux
    ) const;
  void update_flux_speeds();
  void compute_entropy (const Vector<double> &solution,
                        FEValues<dim>        &fe_values,
                        Vector<double>       &entropy) const;
  void compute_entropy_face (const Vector<double> &solution,
                             FEFaceValues<dim>    &fe_values,
                             Vector<double>       &entropy) const;
  void compute_divergence_entropy_flux (const Vector<double> &solution,
                                        FEValues<dim>        &fe_values,
                                        Vector<double>       &divergence) const;
  void compute_velocity(
    const std::vector<double>         & density,
    const std::vector<Tensor<1,dim> > & momentum,
    std::vector<Tensor<1,dim> >       & velocity) const;
  void compute_internal_energy_cell(      std::vector<double>         &internal_energy,
                                    const std::vector<double>         &density,
                                    const std::vector<Tensor<1,dim> > &momentum,
                                    const std::vector<double>         &energy) const;
  void compute_internal_energy_face(      std::vector<double>         &internal_energy,
                                    const std::vector<double>         &density,
                                    const std::vector<Tensor<1,dim> > &momentum,
                                    const std::vector<double>         &energy) const;
  void compute_temperature_cell(      std::vector<double> &temperature,
                                const std::vector<double> &internal_energy) const;
  void compute_temperature_face(      std::vector<double> &temperature,
                                const std::vector<double> &internal_energy) const;
  void compute_pressure(
    const std::vector<double>          & density,
    const std::vector<Tensor<1, dim> > & momentum,
    const std::vector<double>          & energy,
    std::vector<double>                & pressure) const;
  void compute_pressure_cell(      std::vector<double> &pressure,
                                   std::vector<double> &dpdrho,
                                   std::vector<double> &dpdmx,
                                   std::vector<double> &dpdE,
                             const std::vector<double> &density,
                             const std::vector<Tensor<1,dim> > &momentum,
                             const std::vector<double> &energy) const;
  void compute_pressure_face(      std::vector<double> &pressure,
                                   std::vector<double> &dpdrho,
                                   std::vector<double> &dpdmx,
                                   std::vector<double> &dpdE,
                             const std::vector<double> &density,
                             const std::vector<Tensor<1,dim> > &momentum,
                             const std::vector<double> &energy) const;
  void compute_speed_of_sound(      std::vector<double> &speed_of_sound,
                              const std::vector<double> &density,
                              const std::vector<double> &pressure) const;

  double gamma; // gas constant
};

#include "Euler.cc"

#endif
