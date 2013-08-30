/** \file EulerParameters.cc
 *  \brief Provides function definitions for the EulerParameters class.
 */
using namespace dealii;

/** \fn EulerParameters<dim>::EulerParameters()
 *  \brief Constructor for the EulerParameters class.
 */
template<int dim>
EulerParameters<dim>::EulerParameters()
{}

/**
 * \fn    EulerParameters<dim>::declare_euler_parameters(ParameterHandler &parameter_handler)
 * \brief defines input parameters
 * \param parameter_handler parameter handler for the Euler class
 */
template<int dim>
void EulerParameters<dim>::declare_euler_parameters(
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
 * \fn    EulerParameters<dim>::get_euler_parameters(ParameterHandler &parameter_handler)
 * \brief get input parameters from parameter handler
 * \param parameter_handler parameter handler for the Euler class
 */
template<int dim>
void EulerParameters<dim>::get_euler_parameters(
  ParameterHandler &parameter_handler)

{
   // get conservation law parameters
   this->get_conservation_law_parameters(parameter_handler);
   this->n_components = n_euler_components;

   // problem
   parameter_handler.enter_subsection("problem");
   {
      problem_id = parameter_handler.get_integer("problem id");
   }
   parameter_handler.leave_subsection();
}
