/*
 EulerEquations class:
 modified from deal.II step-33 tutorial program
 
 base class for EulerEquationsSinglePhase and EulerEquationsTwoPhase
 contains physical data and flux definitions common to both formulations
*/

#ifndef EulerEquations_h
#define EulerEquations_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/vector.h>

#include "ConservationLaw.h"
#include "EulerEquationsParameters.h"

template <int dim>
class EulerEquations : public ConservationLaw<dim>
{
  public:
    EulerEquations(ParameterHandler &prm);//const std::string &input_filename);

  private:
    // Euler equations parameters
    EulerEquationsParameters<dim> euler_parameters;

    // number of components and position of components in solution vector
    static const unsigned int n_euler_components       = dim + 2;
    static const unsigned int first_momentum_component = 0;
    static const unsigned int density_component        = dim;
    static const unsigned int energy_component         = dim+1;

    // vector of names of each component
    std::vector<std::string>
    get_component_names ();

    // data component interpretation (scalar or vector component) for outputting solution
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_component_interpretations ();

/*
    // compute kinetic energy
    double
    compute_kinetic_energy (const Vector<double> &W);

    // compute pressure
    double
    compute_pressure (const Vector<double> &W);
*/

    void compute_ss_residual(double t, Vector<double> &solution);

    // boundary condition indicators
    enum BoundaryKind
    {
          inflow_boundary,
          outflow_boundary,
          no_penetration_boundary,
          pressure_boundary
    };

/*
    // compute refinement indicators
    void
    compute_refinement_indicators (const DoFHandler<dim> &dof_handler,
                                   const Mapping<dim>    &mapping,
                                   const Vector<double>  &solution,
                                   Vector<double>        &refinement_indicators);
*/

};

#include "EulerEquations.cc"

#endif
