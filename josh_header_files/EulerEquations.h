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
void test_run();

  private:
    // Euler equations parameters
    EulerEquationsParameters<dim> euler_parameters;

    // number of components and position of components in solution vector
    static const unsigned int n_euler_components       = dim + 2;
    static const unsigned int first_momentum_component = 0;
    static const unsigned int density_component        = dim;
    static const unsigned int energy_component         = dim+1;

    // vector of names of each component
//    static
    std::vector<std::string>
    get_component_names ();

    // data component interpretation (scalar or vector component) for outputting solution
//    static
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation ();

    // compute kinetic energy
//    static
    double
    compute_kinetic_energy (const Vector<double> &W);

    // compute pressure
//    static
    double
    compute_pressure (const Vector<double> &W);

    // compute flux matrix f(c)
//    static
    void compute_flux_matrix (const Vector<double> &W,
                              double (&flux)[n_euler_components][dim]);

    // computes numerical normal flux
/*
    template <typename Vector<double>>
    static
    void numerical_normal_flux (const Point<dim>          &normal,
                                const Vector<double>         &Wplus,
                                const Vector<double>         &Wminus,
                                const double               alpha,
                                Sacado::Fad::DFad<double> (&normal_flux)[n_euler_components]);
*/

    // computes forcing vector functions g(c)
//    static
    void compute_forcing_vector (const Vector<double> &W,
                                 double (&forcing)[n_euler_components]);

    // boundary condition indicators
    enum BoundaryKind
    {
          inflow_boundary,
          outflow_boundary,
          no_penetration_boundary,
          pressure_boundary
    };

    // compute refinement indicators
//    static
    void
    compute_refinement_indicators (const DoFHandler<dim> &dof_handler,
                                   const Mapping<dim>    &mapping,
                                   const Vector<double>  &solution,
                                   Vector<double>        &refinement_indicators);

 /*
    class Postprocessor : public DataPostprocessor<dim>
    {
      public:
        Postprocessor (const bool do_schlieren_plot);

        virtual
        void
        compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                           const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                           const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                           const std::vector<Point<dim> >                  &normals,
                                           const std::vector<Point<dim> >                  &evaluation_points,
                                           std::vector<Vector<double> >                    &computed_quantities) const;

        virtual std::vector<std::string> get_names () const;

        virtual
        std::vector<DataComponentInterpretation::DataComponentInterpretation>
        get_data_component_interpretation () const;

        virtual UpdateFlags get_needed_update_flags () const;

      private:
        const bool do_schlieren_plot;
    };
*/
};

#include "EulerEquations.cc"

#endif
