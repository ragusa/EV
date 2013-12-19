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
   velocity_extractor(0)
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
   fe_values[velocity_extractor].get_function_values   (this->current_solution,solution_values);
   fe_values[velocity_extractor].get_function_gradients(this->current_solution,solution_gradients);

   // compute derivative of flux
   std::vector<Tensor<1, dim> > dfdu(this->n_q_points_cell);
   // loop over quadrature points
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      for (int d = 0; d < dim; ++d)
         dfdu[q][d] = solution_values[q];
   
   // NEW VISCOSITY DEFINITION:
   // ===========================================
   /*
   double nu_K = 0.0; // artificial viscosity in cell
   double numerator;  // numerator of artificial viscosity
   // compute cell volume
   double cell_volume = 0.0;
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      cell_volume += fe_values.JxW(q);
   // compute local bilinear form for artificial viscosity
   double b_K = 1.0/(this->dofs_per_cell - 1.0)*cell_volume; // local bilinear form for i != j
   
   // loop over test functions
   for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      // loop over interpolation functions
      for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
         if (j != i)
         {
            numerator = 0.0;
            // loop over quadrature points
            for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
            {
               numerator += dfdu[q]*fe_values[velocity_extractor].gradient(i,q)*
                  fe_values[velocity_extractor].value(j,q)*fe_values.JxW(q);
            }
            nu_K = std::max(nu_K, std::abs(numerator)/(-1.0*b_K));
         }
   std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);
*/
   // ===========================================

   // loop over test functions
   for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
   {
   // ===========================================
   /*
      // compute local bilinear form b_K(u, phi_i) for artificial viscosity
      double b_local = 0.0;
      for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
         if (j != i)
            b_local += this->current_solution(local_dof_indices[j])*b_K;
      cell_residual(i) -= nu_K * b_local;
*/
   // ===========================================
      // loop over quadrature points
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
         cell_residual(i) += (
                               -fe_values[velocity_extractor].value(i,q)
                                *dfdu[q]
                                *solution_gradients[q]
                               -fe_values[velocity_extractor].gradient(i,q)
                                *this->viscosity_cell_q[cell](q)
                                *solution_gradients[q]
                             ) * fe_values.JxW(q);
   }
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
         fe_face_values[velocity_extractor].get_function_gradients   (this->current_solution,solution_gradients_face);

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
            case BurgersParameters<dim>::first_order_1:
            {
               // get max velocity on cell
               std::vector<double> local_solution(this->n_q_points_face);
               fe_face_values.get_function_values(this->current_solution, local_solution);
               double max_velocity = 0.0;
               for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                  max_velocity = std::max( max_velocity, local_solution[q]);

               // compute first-order viscosity
               double cell_diameter = cell->diameter();
               double viscosity_value = burgers_parameters.first_order_viscosity_coef * cell_diameter * max_velocity;
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
               cell_residual(i) += fe_face_values[velocity_extractor].value(i,q)
                                   *viscosity_face(q)
                                   *solution_gradients_face[q]
                                   *fe_face_values.normal_vector(q)
                                   *fe_face_values.JxW(q);
      }
   }
}

/** \fn void Burgers<dim>::compute_ss_Jacobian()
 *  \brief Computes the steady-state Jacobian matrix and stores in system_matrix.
 */
template <int dim>
void Burgers<dim>::compute_ss_Jacobian()
{
   // reset steady-state Jacobian to zero
   this->system_matrix = 0.0;

   // local DoF indices
   std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);

   // velocity at each quadrature point in cell
   std::vector<double> velocity (this->n_q_points_cell);

   // cell matrix
   FullMatrix<double> cell_matrix(this->dofs_per_cell,this->dofs_per_cell);

   // FE values
   FEValues<dim> fe_values(this->fe,this->cell_quadrature,
      update_values | update_gradients | update_JxW_values);

   // ones vector
   Tensor<1,dim> ones_vector;
   for (unsigned int d = 0; d < dim; ++d)
      ones_vector[d] = 1.0;
 
   typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
                                                  endc = this->dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      // reset cell matrix to zero
      cell_matrix = 0;

      cell->get_dof_indices(local_dof_indices);

      fe_values.reinit(cell);
      // get velocity values at each quadrature point in cell
      fe_values[velocity_extractor].get_function_values (this->current_solution, velocity);

      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
         for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
         {
            for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
            {
               cell_matrix(i,j) += fe_values[velocity_extractor].gradient(i,q) * ones_vector
                  * velocity[q] * fe_values[velocity_extractor].value(j,q) * fe_values.JxW(q);
            }
         }
      // add to global matrix
      this->constraints.distribute_local_to_global(cell_matrix, local_dof_indices, this->system_matrix);
   }

}

/** \fn void Burgers<dim>::update_flux_speeds()
 *  \brief Computes the flux speed at each quadrature point in domain and
 *     finds the max in each cell and the max in the entire domain.
 */
template <int dim>
void Burgers<dim>::update_flux_speeds()
{
   FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);
   Tensor<1,dim> dfdu;
   std::vector<double> velocity(this->n_q_points_cell);

   // reset max flux speed
   this->max_flux_speed = 0.0;

   // loop over cells to compute first order viscosity at each quadrature point
   typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
                                                  endc = this->dof_handler.end();
   for (; cell != endc; ++cell)
   {
      fe_values.reinit(cell);
      fe_values[velocity_extractor].get_function_values(this->current_solution, velocity);

      this->max_flux_speed_cell[cell] = 0.0;
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      {
         for (unsigned int d = 0; d < dim; ++d)
            dfdu[d] = velocity[q];
         this->max_flux_speed_cell[cell] = std::max( this->max_flux_speed_cell[cell],
            dfdu.norm() );
      }

      // get max flux speed
      this->max_flux_speed = std::max(this->max_flux_speed, this->max_flux_speed_cell[cell]);
   }
}

/** \fn void Burgers<dim>::compute_entropy(const Vector<double> &solution,
 *                                         FEValues<dim> &fe_values,
 *                                         Vector<double> &entropy) const
 *  \brief Computes entropy at each quadrature point in cell
 *  \param solution solution
 *  \param fe_values FEValues object
 *  \param entropy entropy values at each quadrature point in cell
 */
template <int dim>
void Burgers<dim>::compute_entropy(const Vector<double> &solution,
                                   FEValues<dim>        &fe_values,
                                   Vector<double>       &entropy) const
{
   std::vector<double> velocity(this->n_q_points_cell);
   fe_values[velocity_extractor].get_function_values(solution, velocity);

   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      entropy(q) = 0.5*velocity[q]*velocity[q];
}

/** \fn void Burgers<dim>::compute_entropy_face(const Vector<double> &solution,
 *                                              FEValues<dim> &fe_values_face,
 *                                              Vector<double> &entropy) const
 *  \brief Computes entropy at each quadrature point on face
 *  \param solution solution
 *  \param fe_values_face FEFaceValues object
 *  \param entropy entropy values at each quadrature point on face
 */
template <int dim>
void Burgers<dim>::compute_entropy_face(const Vector<double> &solution,
                                        FEFaceValues<dim>    &fe_values_face,
                                        Vector<double>       &entropy) const
{
   std::vector<double> velocity(this->n_q_points_face);
   fe_values_face[velocity_extractor].get_function_values(solution, velocity);

   for (unsigned int q = 0; q < this->n_q_points_face; ++q)
      entropy(q) = 0.5*velocity[q]*velocity[q];
}

/** \fn void Burgers<dim>::compute_divergence_entropy_flux (const Vector<double> &solution,
 *                                                          FEValues<dim> &fe_values,
 *                                                          Vector<double> &entropy_derivative) const
 *  \brief Computes divergence of entropy flux at each quadrature point in cell
 *  \param solution solution
 *  \param fe_values FEValues object
 *  \param divergence_entropy_flux divergence of entropy flux at each quadrature point in cell
 */
template <int dim>
void Burgers<dim>::compute_divergence_entropy_flux (const Vector<double> &solution,
                                                    FEValues<dim>        &fe_values,
                                                    Vector<double>       &divergence_entropy_flux) const
{
   std::vector<double>         velocity         (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > velocity_gradient(this->n_q_points_cell);

   fe_values[velocity_extractor].get_function_values   (solution, velocity);
   fe_values[velocity_extractor].get_function_gradients(solution, velocity_gradient);

   // constant field v = (1,1,1) (3-D)
   Tensor<1,dim> v;
   for (unsigned int d = 0; d < dim; ++d)
      v[d] = 1.0;

   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
   {
      // compute dot product of constant field v with gradient of u
      double v_dot_velocity_gradient = v * velocity_gradient[q];

      divergence_entropy_flux(q) = velocity[q]*velocity[q]*v_dot_velocity_gradient;
   }
}

/** \fn void Burgers<dim>::output_solution() const
 *  \brief Outputs the solution to .vtk.
 */
template <int dim>
void Burgers<dim>::output_solution () const
{
   if (this->in_final_cycle)
   {
      DataOut<dim> data_out;
      data_out.attach_dof_handler (this->dof_handler);
   
      data_out.add_data_vector (this->current_solution,
                                this->component_names,
                                DataOut<dim>::type_dof_data,
                                this->component_interpretations);
   
      data_out.add_data_vector (this->current_solution, "current_solution");
   
      data_out.build_patches ();
   
      static unsigned int output_file_number = 0;
      if (dim == 1)
      {
         std::string filename = "output/solution-" +
                                Utilities::int_to_string (output_file_number, 3) +
                                ".gpl";
         std::ofstream output (filename.c_str());
         data_out.write_gnuplot (output);
      }
      else
      {
         std::string filename = "output/solution-" +
                                Utilities::int_to_string (output_file_number, 3) +
                                ".vtk";
         std::ofstream output (filename.c_str());
         data_out.write_vtk (output);
      }
   
      ++output_file_number;
   }
}
