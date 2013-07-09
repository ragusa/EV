/*
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/vector.h>

#include "ConservationLaw.h"
#include "EulerEquationsBaseParameters.h"
*/

template <int dim>
EulerEquationsBase<dim>::EulerEquationsBase(const EulerEquationsBaseParameters<dim> &euler_parameters):
   ConservationLaw<dim>(n_euler_components),
   euler_parameters(euler_parameters)
{}

// vector of names of each component
template <int dim>
//    static
std::vector<std::string> EulerEquationsBase<dim>::get_component_names ()
{
   std::vector<std::string> names (dim, "momentum");
   names.push_back ("density");
   names.push_back ("energy_density");

   return names;
}

// data component interpretation (scalar or vector component) for outputting solution
template <int dim>
//    static
std::vector<DataComponentInterpretation::DataComponentInterpretation>
   EulerEquationsBase<dim>::component_interpretation ()
{
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (dim, DataComponentInterpretation::component_is_part_of_vector);
   data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_scalar);
   data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_scalar);

   return data_component_interpretation;
} 
 
// compute kinetic energy
template <int dim>
//    static
double EulerEquationsBase<dim>::compute_kinetic_energy (const Vector<double> &W)
{
   double kinetic_energy = 0.0;
   return kinetic_energy;
}

// compute pressure
template <int dim>
//    static
double EulerEquationsBase<dim>::compute_pressure (const Vector<double> &W)
{
   double pressure = 0.0;
   return pressure;
}

// compute flux matrix f(c)
template <int dim>
//    static
void EulerEquationsBase<dim>::compute_flux_matrix (const Vector<double> &W,
                                                   double (&flux)[n_euler_components][dim])
{}

/*
// computes numerical normal flux
template <int dim>
    template <typename InputVector>
    static
    void numerical_normal_flux (const Point<dim>          &normal,
                                const InputVector         &Wplus,
                                const InputVector         &Wminus,
                                const double               alpha,
                                Sacado::Fad::DFad<double> (&normal_flux)[n_euler_components]);
*/

// computes forcing vector functions g(c)
template <int dim>
//    static
void EulerEquationsBase<dim>::compute_forcing_vector (const Vector<double> &W,
                                                      double (&forcing)[n_euler_components])
{}

// compute refinement indicators
template <int dim>
//    static
void EulerEquationsBase<dim>::compute_refinement_indicators (const DoFHandler<dim> &dof_handler,
                                                             const Mapping<dim>    &mapping,
                                                             const Vector<double>  &solution,
                                                             Vector<double>        &refinement_indicators)
{}
