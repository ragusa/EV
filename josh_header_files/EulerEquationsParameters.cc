/** \file EulerEquationsParameters.cc
 *  \brief Provides the function definitions for the EulerEquationsParameters class.
 */
using namespace dealii;

/** \fn EulerEquationsParameters<dim>::EulerEquationsParameters()
 *  \brief Constructor for the EulerEquationsParametrs class.
 */
template<int dim>
EulerEquationsParameters<dim>::EulerEquationsParameters():
   initial_conditions_expressions(n_components,"0")
{}

/**
 * \fn    EulerEquationsParameters<dim>::declare_parameters(ParameterHandler &parameter_handler)
 * \brief Declares input parameters for the EulerEquations class.
 * \param parameter_handler parameter handler for EulerEquations class
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
 * \fn    EulerEquationsParameters<dim>::get_parameters(ParameterHandler &parameter_handler)
 * \brief get input parameters from parameter handler
 * \param parameter_handler parameter handler for EulerEquations class
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
      //std::vector<std::string> expressions (n_components,"0.0");
      for (int c = 0; c < n_components; c++)
          initial_conditions_expressions[c] = parameter_handler.get("initial conditions " + Utilities::int_to_string(c));
   }
   parameter_handler.leave_subsection();
}
