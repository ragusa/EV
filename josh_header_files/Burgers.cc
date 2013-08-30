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
{} 

/** \fn std::vector<std::string> Burgers<dim>::get_component_names()
 *  \brief Returns the names of each component.
 *  \return vector of names of each component
 */
template <int dim>
std::vector<std::string> Burgers<dim>::get_component_names ()
{
   std::vector<std::string> names(1,"velocity");
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

/** \fn void Burgers<dim>::define_problem()
 *  \brief Create the domain, compute volume, define initial
 *         conditions, and define boundary conditions and exact
 *         solution if it exists.
 */
template <int dim>
void Burgers<dim>::define_problem()
{
   switch (burgers_parameters.problem_id)
   {
      case 0: // 1-D, Dirichlet boundary conditions, sin(2*pi*x)
      {
         Assert(dim==1,ExcImpossibleInDim(dim));
         double domain_start = 0;
         double domain_width = 1.0;
         this->domain_volume = std::pow(domain_width,dim);
         GridGenerator::hyper_cube(this->triangulation, domain_start, domain_start + domain_width);
         // only 1 type of BC: zero Dirichlet; leave boundary indicators as zero
         this->n_boundaries = 1;
         // set all boundary indicators to zero
         typename Triangulation<dim>::cell_iterator cell = this->triangulation.begin(),
                                                    endc = this->triangulation.end();
         for (; cell != endc; ++cell)
            for (unsigned int face = 0; face < this->faces_per_cell; ++face)
               if (cell->face(face)->at_boundary())
                  cell->face(face)->set_boundary_indicator(0);

         this->boundary_types.resize(this->n_boundaries);
         this->boundary_types[0].resize(this->n_components);
         this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet;
         this->dirichlet_function_strings.resize(this->n_components);
         this->dirichlet_function_strings[0] = "0";
         this->use_exact_solution_as_BC = false;
         // initial conditions
         this->initial_conditions_strings[0] = "sin(2*pi*x)";
         // exact solution
         this->has_exact_solution = false;
         break;
      }
      case 1: // Guermond 2-d test problem
      {
         Assert(dim==2,ExcImpossibleInDim(dim));
         double domain_start = 0;
         double domain_width = 1.0;
         this->domain_volume = std::pow(domain_width,dim);
         GridGenerator::hyper_cube(this->triangulation, domain_start, domain_start + domain_width);
         // only 1 type of BC: Dirichlet with exact solution
         this->n_boundaries = 1;
         // set all boundary indicators to zero
         typename Triangulation<dim>::cell_iterator cell = this->triangulation.begin(),
                                                    endc = this->triangulation.end();
         for (; cell != endc; ++cell)
            for (unsigned int face = 0; face < this->faces_per_cell; ++face)
               if (cell->face(face)->at_boundary())
                  cell->face(face)->set_boundary_indicator(0);

         this->boundary_types.resize(this->n_boundaries);
         this->boundary_types[0].resize(this->n_components);
         this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet;
         this->use_exact_solution_as_BC = true;
         // initial conditions
         this->initial_conditions_strings[0] =  "if(x<0.5,";
         this->initial_conditions_strings[0] +=    "if(y>0.5,";
         this->initial_conditions_strings[0] +=       "-0.2,0.5),";
         this->initial_conditions_strings[0] +=    "if(y<0.5,";
         this->initial_conditions_strings[0] +=       "0.8,-1))";
         // exact solution
         this->has_exact_solution = true;
         this->exact_solution_strings[0] =  "if(x<0.5-0.6*t,";
         this->exact_solution_strings[0] +=    "if(y>0.5+0.15*t,";
         this->exact_solution_strings[0] +=       "-0.2,0.5),";
         this->exact_solution_strings[0] +=  "if(x<0.5-0.25*t,";
         this->exact_solution_strings[0] +=    "if(y>-8./7.*x+15./14.-15./28.*t,";
         this->exact_solution_strings[0] +=       "-1.0,0.5),";
         this->exact_solution_strings[0] +=  "if(x<0.5+0.5*t,";
         this->exact_solution_strings[0] +=    "if(y>x/6.+5./12.-5./24.*t,";
         this->exact_solution_strings[0] +=       "-1.0,0.5),";
         this->exact_solution_strings[0] +=  "if(x<0.5+0.8*t,";
         this->exact_solution_strings[0] +=    "if(y>x-5./(18.*t)*(x+t-0.5)^2,";
         this->exact_solution_strings[0] +=       "-1.0,(2*x-1)/(2.*t)),";
         this->exact_solution_strings[0] +=  "if(y>0.5-0.1*t,";
         this->exact_solution_strings[0] +=       "-1.0,0.8)";
         this->exact_solution_strings[0] +=  "))))";
         break;
      }
      default:
         Assert(false,ExcNotImplemented());
   }
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
               for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                  viscosity_face(q) = burgers_parameters.constant_viscosity_value;
               break;
            }
            case BurgersParameters<dim>::first_order:
            {
               // get max velocity on cell
               std::vector<double> local_solution(this->n_q_points_face);
               fe_face_values.get_function_values(this->current_solution, local_solution);
               double max_velocity = 0.0;
               for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                  max_velocity = std::max( max_velocity, local_solution[q]);

               // compute first-order viscosity
               double cell_diameter = cell->diameter();
               double viscosity_value = 0.5 * burgers_parameters.first_order_viscosity_coef * cell_diameter * max_velocity;
               for (unsigned int q = 0; q < this->n_q_points_face; ++q)
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
         for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
            // loop over quadrature points
            for (unsigned int q = 0; q < this->n_q_points_face; ++q)
               cell_residual(i) += fe_face_values[velocity].value(i,q)
                                   *viscosity_face(q)
                                   *solution_gradients_face[q]
                                   *fe_face_values.normal_vector(q)
                                   *fe_face_values.JxW(q);
      }
   }
}

/** \fn Tensor<1,dim> Burgers<dim>::flux_derivative(const double u)
 *  \brief Computes the derivative of the flux function with respect to
 *         the solution vector.
 *  \param u solution at a point
 *  \return derivative of the flux with respect to u
 */
template <int dim>
Tensor<1,dim> Burgers<dim>::flux_derivative(const double u)
{
   Tensor<1,dim> dfdu;
   
   for (unsigned int d = 0; d < dim; ++d)
      dfdu[d] = u;
   
   return dfdu;
}

/** \fn double Burgers<dim>::entropy(const double u) const
 *  \brief Computes entropy at a point.
 *  \param u solution at a point
 *  \return entropy at a point
 */
template <int dim>
double Burgers<dim>::entropy(const double u) const
{
   return 0.5*std::pow(u,2);
}

/** \fn double Burgers<dim>::entropy_derivative(const double u) const
 *  \brief Computes derivative of entropy at a point.
 *  \param u solution at a point
 *  \return derivative of entropy at a point
 */
template <int dim>
double Burgers<dim>::entropy_derivative(const double u) const
{
   return u;
}
