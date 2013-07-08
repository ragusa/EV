/*
 EulerEquationsBase class:
 modified from deal.II step-33 tutorial program
 
 base class for EulerEquationsSinglePhase and EulerEquationsTwoPhase
 contains physical data and flux definitions common to both formulations
*/

#ifndef EulerEquationsBase_h
#define EulerEquationsBase_h

#include "ConservationLaw.h"
#include "EulerEquationsBaseParameters.h"

template <int dim>
class EulerEquationsBase : public ConservationLaw
{
  public:
    EulerEquationsBase(const EulerEquationsBaseParameters &parameters);

  private:
    // number of components and position of components in solution vector
    static const unsigned int n_components             = dim + 2;
    static const unsigned int first_momentum_component = 0;
    static const unsigned int density_component        = dim;
    static const unsigned int energy_component         = dim+1;

    // vector of names of each component
    static
    std::vector<std::string>
    component_names ();

    // data component interpretation (scalar or vector component) for outputting solution
    static
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation ();

    // gas constant
    static const double gas_gamma;

    // compute kinetic energy
    template <typename number, typename InputVector>
    static
    number
    compute_kinetic_energy (const InputVector &W);

    // compute pressure
    template <typename number, typename InputVector>
    static
    number
    compute_pressure (const InputVector &W);

    // compute flux matrix f(c)
    template <typename InputVector, typename number>
    static
    void compute_flux_matrix (const InputVector &W,
                              number (&flux)[n_components][dim]);

    // computes numerical normal flux
    template <typename InputVector>
    static
    void numerical_normal_flux (const Point<dim>          &normal,
                                const InputVector         &Wplus,
                                const InputVector         &Wminus,
                                const double               alpha,
                                Sacado::Fad::DFad<double> (&normal_flux)[n_components]);

    // computes forcing vector functions g(c)
    template <typename InputVector, typename number>
    static
    void compute_forcing_vector (const InputVector &W,
                                 number (&forcing)[n_components]);

    // boundary condition indicators
    enum BoundaryKind
    {
          inflow_boundary,
          outflow_boundary,
          no_penetration_boundary,
          pressure_boundary
    };

    // compute W_minus
    template <typename DataVector>
    static
    void
    compute_Wminus (const BoundaryKind  (&boundary_kind)[n_components],
                    const Point<dim>     &normal_vector,
                    const DataVector     &Wplus,
                    const Vector<double> &boundary_values,
                    const DataVector     &Wminus);

    // compute refinement indicators
    static
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

#include "EulerEquationsBase.cc"

#endif
