/*
template <int dim>
void EulerEquations<dim>::print_test_message_from_derived_class()
{
   this->print_test_message_from_base_class();
}
*/

template <int dim>
EulerEquations<dim>::EulerEquations(ParameterHandler &prm):
   ConservationLaw<dim>(prm,n_euler_components)
{
   // get Euler parameters
   euler_parameters.get_parameters(prm);
   // initialize initial conditions
   this->initial_conditions.initialize (FunctionParser<dim>::default_variable_names(),
                                        euler_parameters.initial_conditions_expressions,
                                        std::map<std::string, double>());
   this->component_names = get_component_names();
   this->component_interpretations = get_component_interpretations();
}

// vector of names of each component
template <int dim>
std::vector<std::string> EulerEquations<dim>::get_component_names ()
{
   std::vector<std::string> names (dim, "momentum");
   names.push_back ("density");
   names.push_back ("energy_density");

   return names;
}

// data component interpretation (scalar or vector component) for outputting solution
template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
   EulerEquations<dim>::get_component_interpretations ()
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
double EulerEquations<dim>::compute_kinetic_energy (const Vector<double> &W)
{
   double kinetic_energy = 0.0;
   return kinetic_energy;
}

// compute pressure
template <int dim>
//    static
double EulerEquations<dim>::compute_pressure (const Vector<double> &W)
{
   double pressure = 0.0;
   return pressure;
}

// compute flux matrix f(c)
template <int dim>
//    static
void EulerEquations<dim>::compute_flux_matrix (const Vector<double> &W,
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
void EulerEquations<dim>::compute_forcing_vector (const Vector<double> &W,
                                                      double (&forcing)[n_euler_components])
{}

// compute refinement indicators
template <int dim>
//    static
void EulerEquations<dim>::compute_refinement_indicators (const DoFHandler<dim> &dof_handler,
                                                             const Mapping<dim>    &mapping,
                                                             const Vector<double>  &solution,
                                                             Vector<double>        &refinement_indicators)
{}
