/** \file BurgersParameters.cc
 *  \brief Provides function definitions for the BurgersParameters class.
 */
using namespace dealii;

/** \fn BurgersParameters<dim>::BurgersParameters()
 *  \brief Constructor for the BurgersParameters class.
 */
template<int dim>
BurgersParameters<dim>::BurgersParameters():
   initial_conditions_expressions(n_burgers_components,"0")
{}

/**
 * \fn    BurgersParameters<dim>::declare_burgers_parameters(ParameterHandler &parameter_handler)
 * \brief defines input parameters
 * \param parameter_handler parameter handler for the Burgers class
 */
template<int dim>
void BurgersParameters<dim>::declare_burgers_parameters(
  ParameterHandler &parameter_handler)
{
   // declare conservation law parameters
   ConservationLawParameters<dim>::declare_conservation_law_parameters(parameter_handler);

   // initial conditions
   parameter_handler.enter_subsection("initial conditions");
   {
      for (int c = 0; c < n_burgers_components; ++c)
         parameter_handler.declare_entry("initial conditions " + Utilities::int_to_string(c),
         "0.0",
         Patterns::Anything(),
         "initial conditions for component computed from x,y,z");
   }
   parameter_handler.leave_subsection();
}

/**
 * \fn    BurgersParameters<dim>::get_burgers_parameters(ParameterHandler &parameter_handler)
 * \brief get input parameters from parameter handler
 * \param parameter_handler parameter handler for the Burgers class
 */
template<int dim>
void BurgersParameters<dim>::get_burgers_parameters(
  ParameterHandler &parameter_handler)

{
   // get conservation law parameters
   this->get_conservation_law_parameters(parameter_handler);
   this->n_components = n_burgers_components;

   // initial conditions
   parameter_handler.enter_subsection("initial conditions");
   {
      //std::vector<std::string> expressions (n_components,"0.0");
      for (int c = 0; c < n_burgers_components; c++)
          initial_conditions_expressions[c] = parameter_handler.get("initial conditions " + Utilities::int_to_string(c));
   }
   parameter_handler.leave_subsection();
}
