using namespace dealii;

template<int dim>
EulerEquationsParameters<dim>::EulerEquationsParameters():
   initial_conditions(n_components)
{}

/**
 * \fn    EulerEquationsParameters::declare_parameters
 * \brief defines input parameters
 */
template<int dim>
void EulerEquationsParameters<dim>::declare_parameters(
  ParameterHandler &parameter_handler)
{
   parameter_handler.declare_entry("First input parameter", "5.0",
      Patterns::Double(),"Description of first input parameter");
   parameter_handler.declare_entry("Second input parameter", "5.0",
      Patterns::Double(),"Description of second input parameter");
   parameter_handler.declare_entry("Third input parameter", "5.0",
      Patterns::Double(),"Description of third input parameter");

   // initial conditions
   parameter_handler.enter_subsection("initial conditions");
   {
      for (int c = 0; c < n_components; ++c)
         parameter_handler.declare_entry("initial conditions " + Utilities::int_to_string(c),
         "0.0",
         Patterns::Anything(),
         "initial conditions for component computed from x,y,z");
   }
   parameter_handler.leave_subsection();
}

/**
 * \fn    EulerEquationsParameters::get_parameters
 * \brief get input parameters from parameter handler
 */
template<int dim>
void EulerEquationsParameters<dim>::get_parameters(
  ParameterHandler &parameter_handler)
{
   input1 = parameter_handler.get_double("First input parameter");
   input2 = parameter_handler.get_double("Second input parameter");
   input3 = parameter_handler.get_double("Third input parameter");

   // initial conditions
   parameter_handler.enter_subsection("initial conditions");
   {
      std::vector<std::string> expressions (n_components,"0.0");
      for (int c = 0; c < n_components; c++)
          expressions[c] = parameter_handler.get("initial conditions " + Utilities::int_to_string(c));
      initial_conditions.initialize (FunctionParser<dim>::default_variable_names(),
                                     expressions,
                                     std::map<std::string, double>());
   }
   parameter_handler.leave_subsection();
}
