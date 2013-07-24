using namespace dealii;

template<int dim>
BurgersParameters<dim>::BurgersParameters():
   initial_conditions_expressions(n_components,"0")
{}

/**
 * \fn    BurgersParameters::declare_parameters
 * \brief defines input parameters
 */
template<int dim>
void BurgersParameters<dim>::declare_parameters(
  ParameterHandler &parameter_handler)
{
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
 * \fn    BurgersParameters::get_parameters
 * \brief get input parameters from parameter handler
 */
template<int dim>
void BurgersParameters<dim>::get_parameters(
  ParameterHandler &parameter_handler)
{
   // initial conditions
   parameter_handler.enter_subsection("initial conditions");
   {
      //std::vector<std::string> expressions (n_components,"0.0");
      for (int c = 0; c < n_components; c++)
          initial_conditions_expressions[c] = parameter_handler.get("initial conditions " + Utilities::int_to_string(c));
   }
   parameter_handler.leave_subsection();
}
