/** \file BurgersParameters.cc
 *  \brief Provides function definitions for the BurgersParameters class.
 */
using namespace dealii;

/** \fn BurgersParameters<dim>::BurgersParameters()
 *  \brief Constructor for the BurgersParameters class.
 */
template<int dim>
BurgersParameters<dim>::BurgersParameters()
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

   // problem
   parameter_handler.enter_subsection("problem");
   {
      parameter_handler.declare_entry("problem id",
                                      "0",
                                      Patterns::Integer(),
                                      "ID for description of the problem");
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

   // problem
   parameter_handler.enter_subsection("problem");
   {
      problem_id = parameter_handler.get_integer("problem id");
   }
   parameter_handler.leave_subsection();
}
