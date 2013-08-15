/** \file Burgers.cc
 *  \brief Provides the function definitions for the Burgers class.
 */

/** \fn Burgers<dim>::Burgers(const BurgersParameters<dim> &params)
 *  \brief Constructor for the Burgers class.
 *  \param params Burgers' equation parameters
 */
template <int dim>
Burgers<dim>::Burgers(const BurgersParameters<dim> &params):
   ConservationLaw<dim>(params),
   burgers_parameters(params),
   velocity(0)
{
   // initialize initial conditions
   std::map<std::string,double> constants;
   constants["pi"] = numbers::PI;

   this->initial_conditions.initialize (FunctionParser<dim>::default_variable_names(),
                                        burgers_parameters.initial_conditions_expressions,
                                        constants);
   this->component_names           = this->get_component_names();
   this->component_interpretations = this->get_component_interpretations();
}

/** \fn std::vector<std::string> Burgers<dim>::get_component_names()
 *  \brief Returns the names of each component.
 *  \return vector of names of each component
 */
template <int dim>
std::vector<std::string> Burgers<dim>::get_component_names ()
{
   std::vector<std::string> names (1, "velocity");

   return names;
}

/** \fn std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *      Burgers<dim>::get_component_interpretations()
 *  \brief Returns the interpretations for each component.
 *
 *  This function returns the interpretation of each component,
 *  i.e., whether each component is a scalar or a component of
 *  a vector.
 *  \return data component interpretations
 */
template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
   Burgers<dim>::get_component_interpretations ()
{
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (1, DataComponentInterpretation::component_is_scalar);

   return data_component_interpretation;
} 

/** \fn void Burgers<dim>::compute_cell_ss_residual()
 *  \brief Computes the contribution of the steady-state residual
 *         from the current cell.
 */
template <int dim>
void Burgers<dim>::compute_cell_ss_residual(FEValues<dim> &fe_values,
                                            const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            Vector<double> &cell_residual)
{
   // reinitialize fe values for cell
   fe_values.reinit(cell);

   // get current solution values and gradients
   std::vector<double>          solution_values   (this->n_q_points_cell);
   std::vector<Tensor<1, dim> > solution_gradients(this->n_q_points_cell);
   fe_values[velocity].get_function_values   (this->current_solution,solution_values);
   fe_values[velocity].get_function_gradients(this->current_solution,solution_gradients);

   // compute derivative of flux
   std::vector<Tensor<1, dim> > dfdu(this->n_q_points_cell);
   // loop over quadrature points
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      dfdu[q] = flux_derivative(solution_values[q]);
   
   // loop over test functions
   for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      // loop over quadrature points
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
         cell_residual(i) -= (
                                fe_values[velocity].value(i,q)
                                *dfdu[q]
                                *solution_gradients[q]
                             +  fe_values[velocity].gradient(i,q)
                                *this->viscosity_cell_q[cell](q)
                                *solution_gradients[q]
                             ) * fe_values.JxW(q);
}

/** \fn void Burgers<dim>::compute_face_ss_residual()
 *  \brief Computes the contribution of the steady-state residual
 *         from the faces of the current cell.
 */
template <int dim>
void Burgers<dim>::compute_face_ss_residual(FEFaceValues<dim> &fe_face_values,
                                            const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            Vector<double> &cell_residual)
{
/*
   // loop over faces
   for (unsigned int face = 0; face < this->faces_per_cell; ++face)
   {
      // add term for boundary faces
      if (cell->at_boundary(face))
      {
         fe_face_values.reinit(cell, face);

         std::vector<Tensor<1, dim> > solution_gradients_face(this->n_q_points_face);
         fe_face_values[velocity].get_function_gradients   (this->current_solution,solution_gradients_face);

         // compute viscosity
         Vector<double> viscosity_face(this->n_q_points_face);
         switch (burgers_parameters.viscosity_type)
         {
            case BurgersParameters<dim>::constant:
            {
               for (unsigned int q = 0; q < n_q_points_face; ++q)
                  viscosity_face(q) = burgers_parameters.constant_viscosity_value;
               break;
            }
            case BurgersParameters<dim>::first_order:
            {
               // get max velocity on cell
               std::vector<double> local_solution(n_q_points_face);
               fe_face_values.get_function_values(this->current_solution, local_solution);
               double max_velocity = 0.0;
               for (unsigned int q = 0; q < n_q_points_face; ++q)
                  max_velocity = std::max( max_velocity, local_solution[q]);

               // compute first-order viscosity
               double cell_diameter = cell->diameter();
               double viscosity_value = 0.5 * burgers_parameters.first_order_viscosity_coef * cell_diameter * max_velocity;
               for (unsigned int q = 0; q < n_q_points_face; ++q)
                  viscosity_face(q) = viscosity_value;
               
               break;
            }
            default:
            {
               Assert(false,ExcNotImplemented());
               break;
            }
         }
         // loop over test functions
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
            // loop over quadrature points
            for (unsigned int q = 0; q < n_q_points_face; ++q)
               cell_residual(i) += fe_face_values[velocity].value(i,q)
                                   *viscosity_face(q)
                                   *solution_gradients_face[q]
                                   *fe_face_values.normal_vector(q)
                                   *fe_face_values.JxW(q);
      }
   }
*/
}

template <int dim>
Tensor<1,dim> Burgers<dim>::flux_derivative(const double u)
{
   Tensor<1,dim> dfdu;
   
   for (unsigned int d = 0; d < dim; ++d)
      dfdu[d] = u;
   
   return dfdu;
}

template <int dim>
double Burgers<dim>::entropy(const double u) const
{
   return 0.5*std::pow(u,2);
}

template <int dim>
double Burgers<dim>::entropy_derivative(const double u) const
{
   return u;
}
