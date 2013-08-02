/** \file BurgersParameters.cc
 *  \brief Provides function definitions for the BurgersParameters class.
 */
using namespace dealii;

/** \fn BurgersParameters<dim>::BurgersParameters()
 *  \brief Constructor for the BurgersParameters class.
 */
template<int dim>
BurgersParameters<dim>::BurgersParameters():
   initial_conditions_expressions(n_components,"0")
{}

/**
 * \fn    BurgersParameters<dim>::declare_parameters(ParameterHandler &parameter_handler)
 * \brief defines input parameters
 * \param parameter_handler parameter handler for the Burgers class
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

   // artificial viscosity
   parameter_handler.enter_subsection("artificial viscosity");
   {
      parameter_handler.declare_entry("viscosity type",
                                      "constant",
                                      Patterns::Anything(),
                                      "choice for artificial viscosity");
      parameter_handler.declare_entry("constant viscosity value",
                                      "1e-3",
                                      Patterns::Double(),
                                      "viscosity value if constant viscosity chosen");
      parameter_handler.declare_entry("first order viscosity coefficient",
                                      "1e-3",
                                      Patterns::Double(),
                                      "tuning constant value to be used with first-order viscosity");
      parameter_handler.declare_entry("entropy viscosity coefficient",
                                      "1e-3",
                                      Patterns::Double(),
                                      "tuning constant value to be used with entropy viscosity");
   }
   parameter_handler.leave_subsection();
 
}

/**
 * \fn    BurgersParameters<dim>::get_parameters(ParameterHandler &parameter_handler)
 * \brief get input parameters from parameter handler
 * \param parameter_handler parameter handler for the Burgers class
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

   // artificial viscosity
   parameter_handler.enter_subsection("artificial viscosity");
   {
      const std::string viscosity_choice = parameter_handler.get("viscosity type");
      if (viscosity_choice == "none")
         viscosity_type = none;
      else if (viscosity_choice == "constant")
         viscosity_type = constant;
      else if (viscosity_choice == "first_order")
         viscosity_type = first_order;
      else if (viscosity_choice == "entropy")
         viscosity_type = entropy;
      else
         Assert(false,ExcNotImplemented());

      constant_viscosity_value = parameter_handler.get_double("constant viscosity value");
      first_order_viscosity_coef = parameter_handler.get_double("first order viscosity coefficient");
      entropy_viscosity_coef = parameter_handler.get_double("entropy viscosity coefficient");
   }
   parameter_handler.leave_subsection();
}
