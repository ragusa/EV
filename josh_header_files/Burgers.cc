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
   burgers_parameters(params)
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

/** \fn Burgers<dim>::compute_ss_residual(double t, Vector<double> &solution)
 *  \brief Computes the steady-state residual for Burgers' equation.
 *
 *  This function computes the steady-state residual \f$\mathbf{f_{ss}}\f$ for the conservation law
 *  \f[
 *    \frac{\partial\mathbf{u}}{\partial t} 
 *    + \nabla \cdot \mathbf{f}(\mathbf{u}) = \mathbf{g}(\mathbf{u}),
 *  \f]
 *  which for component \f$i\f$ is
 *  \f[
 *    \mathbf{f_{ss}} = (\mathbf{\psi}, -\nabla \cdot \mathbf{f}(\mathbf{u}) + \mathbf{g}(\mathbf{u}))_\Omega.
 *  \f]
 *  For vicous Burgers' equation, this is the following:
 *  \f[
 *    \mathbf{f_{ss}} = -(\mathbf{\psi},u u_x)_\Omega + (\mathbf{\psi},\nu u_{xx})_\Omega .
 *  \f]
 *  After integration by parts, this is
 *  \f[
 *    \mathbf{f_{ss}} = -(\mathbf{\psi},u u_x)_\Omega - (\mathbf{{\psi}_x},\nu u_{x})_\Omega 
 *    + (\mathbf{\psi},\nu u_{x})_{\partial\Omega}.
 *  \f]
 *  \param t time at which the steady-state residual is to be evaluated
 *  \param solution at which to evaluate the steady-state residual
 */
template <int dim>
void Burgers<dim>::compute_ss_residual(double t, Vector<double> &solution)
{
   const FEValuesExtractors::Scalar velocity (0);

   FEValues<dim> fe_values (this->fe, this->quadrature,
                            update_values | update_gradients | update_JxW_values);
   FEFaceValues<dim> fe_face_values (this->fe, this->face_quadrature,
                            update_values | update_gradients | update_normal_vectors | update_JxW_values);

   const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
   const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;
   const unsigned int n_q_points_cell    = this->quadrature.size();
   const unsigned int n_q_points_face    = this->face_quadrature.size();

   Vector<double> cell_residual(dofs_per_cell);

   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
                                                  endc = this->dof_handler.end();
   // loop over cells
   for (; cell!=endc; ++cell)
   {
      cell_residual = 0;
      fe_values.reinit(cell);
      // get global DoF indices
      cell->get_dof_indices(local_dof_indices);

      std::vector<double>          solution_values   (n_q_points_cell);
      std::vector<Tensor<1, dim> > solution_gradients(n_q_points_cell);
      fe_values[velocity].get_function_values   (solution,solution_values);
      fe_values[velocity].get_function_gradients(solution,solution_gradients);

      std::vector<Tensor<1, dim> > flux_derivative(n_q_points_cell);
      // loop over quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
      {
         Tensor<1, dim> values_q;
         for (int d = 0; d < dim; ++d)
            values_q[d] = solution_values[q];
         flux_derivative[q] = values_q;
      }
      
      // loop over test functions
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
         // loop over quadrature points
         for (unsigned int q = 0; q < n_q_points_cell; ++q)
            cell_residual(i) -= fe_values[velocity].value(i,q)
                                *flux_derivative[q]
                                *solution_gradients[q]
                                *fe_values.JxW(q);

      // add artificial viscosity if requested
      if (burgers_parameters.viscosity_type != BurgersParameters<dim>::none)
      {
         // compute viscosity
         Vector<double> viscosity(n_q_points_cell);
         switch (burgers_parameters.viscosity_type)
         {
            case BurgersParameters<dim>::constant:
            {
               for (unsigned int q = 0; q < n_q_points_cell; ++q)
                  viscosity(q) = burgers_parameters.constant_viscosity_value;
               break;
            }
            case BurgersParameters<dim>::first_order:
            {
               // get max velocity on cell
               std::vector<double> local_solution(n_q_points_cell);
               fe_values.get_function_values(solution, local_solution);
               double max_velocity = 0.0;
               for (unsigned int q = 0; q < n_q_points_cell; ++q)
                  max_velocity = std::max( max_velocity, local_solution[q]);

               // compute first-order viscosity
               double cell_diameter = cell->diameter();
               double viscosity_value = 0.5 * burgers_parameters.first_order_viscosity_coef * cell_diameter * max_velocity;
               for (unsigned int q = 0; q < n_q_points_cell; ++q)
                  viscosity(q) = viscosity_value;
               
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
            for (unsigned int q = 0; q < n_q_points_cell; ++q)
               cell_residual(i) -= fe_values[velocity].gradient(i,q)
                                   *viscosity(q)
                                   *solution_gradients[q]
                                   *fe_values.JxW(q);

         // loop over faces
         for (unsigned int face = 0; face < faces_per_cell; ++face)
         {
            // add term for boundary faces
            if (cell->at_boundary(face))
            {
               fe_face_values.reinit(cell, face);

               std::vector<Tensor<1, dim> > solution_gradients_face(n_q_points_face);
               fe_face_values[velocity].get_function_gradients   (solution,solution_gradients_face);

               // compute viscosity
               Vector<double> viscosity_face(n_q_points_face);
               switch (burgers_parameters.viscosity_type)
               {
                  case BurgersParameters<dim>::constant:
                  {
                     for (unsigned int q = 0; q < n_q_points_face; ++q)
                        viscosity(q) = burgers_parameters.constant_viscosity_value;
                     break;
                  }
                  case BurgersParameters<dim>::first_order:
                  {
                     // get max velocity on cell
                     std::vector<double> local_solution(n_q_points_face);
                     fe_face_values.get_function_values(solution, local_solution);
                     double max_velocity = 0.0;
                     for (unsigned int q = 0; q < n_q_points_face; ++q)
                        max_velocity = std::max( max_velocity, local_solution[q]);
      
                     // compute first-order viscosity
                     double cell_diameter = cell->diameter();
                     double viscosity_value = 0.5 * burgers_parameters.first_order_viscosity_coef * cell_diameter * max_velocity;
                     for (unsigned int q = 0; q < n_q_points_face; ++q)
                        viscosity(q) = viscosity_value;
                     
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
      }

      // aggregate local residual into global residual
      // loop over test functions
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
         this->ss_residual(local_dof_indices[i]) += cell_residual(i);
   }

   // apply boundary conditions: zero Dirichlet
   std::map<unsigned int, double> boundary_values;
   VectorTools::interpolate_boundary_values(this->dof_handler,
                                            0,
                                            ZeroFunction<dim>(),
                                            boundary_values);
   for (std::map<unsigned int, double>::iterator bv = boundary_values.begin(); bv != boundary_values.end(); ++bv)
      this->ss_residual(bv->first) = (bv->second);
}
