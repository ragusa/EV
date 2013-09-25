/** \file Euler.cc
 *  \brief Provides the function definitions for the Euler class.
 */

/** \fn Euler<dim>::Euler(const EulerParameters<dim> &params)
 *  \brief Constructor for the Euler class.
 *  \param params Euler equation parameters
 */
template <int dim>
Euler<dim>::Euler(const EulerParameters<dim> &params):
   ConservationLaw<dim>(params),
   euler_parameters(params),
   density_extractor(0),
   momentum_extractor(1),
   energy_extractor(1+dim)
{} 

/** \fn std::vector<std::string> Euler<dim>::get_component_names()
 *  \brief Returns the names of each component.
 *  \return vector of names of each component
 */
template <int dim>
std::vector<std::string> Euler<dim>::get_component_names ()
{
   std::vector<std::string> names(n_euler_components);

   names[0] = "density";
   for (int d = 0; d < dim; ++d) names[1+d] = "momentum";
   names[dim+1] = "energy";

   return names;
}

/** \fn std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *      Euler<dim>::get_component_interpretations()
 *  \brief Returns the interpretations for each component.
 *
 *  This function returns the interpretation of each component,
 *  i.e., whether each component is a scalar or a component of
 *  a vector.
 *  \return data component interpretations
 */
template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
   Euler<dim>::get_component_interpretations ()
{
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
      component_interpretations(n_euler_components);
   
   component_interpretations[0] = DataComponentInterpretation::component_is_scalar;
   for (int d = 0; d < dim; ++d)
      component_interpretations[1+d] = DataComponentInterpretation::component_is_part_of_vector;
   component_interpretations[dim+1] = DataComponentInterpretation::component_is_scalar;

   return component_interpretations;
} 

/** \fn void Euler<dim>::define_problem()
 *  \brief Create the domain, compute volume, define initial
 *         conditions, and define boundary conditions and exact
 *         solution if it exists.
 */
template <int dim>
void Euler<dim>::define_problem()
{
   switch (euler_parameters.problem_id)
   {
      case 0: // 1-D Shock tube problem
      {
         Assert(dim==1,ExcImpossibleInDim(dim));
         double domain_start = 0;
         double domain_width = 1.0;
         this->domain_volume = std::pow(domain_width,dim);
         GridGenerator::hyper_cube(this->triangulation, domain_start, domain_start + domain_width);
         // only 1 type of BC: zero Dirichlet
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
         this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet; // density has Dirichlet BC
         this->boundary_types[0][1] = ConservationLaw<dim>::dirichlet; // x-momentum has Dirichlet BC
         this->boundary_types[0][2] = ConservationLaw<dim>::dirichlet; // energy has Dirichlet BC
         this->dirichlet_function_strings.resize(this->n_components);
         this->dirichlet_function_strings[0] = "if(x<0.5,1.0,0.125)"; // BC for density
         this->dirichlet_function_strings[1] = "0";                   // BC for x-momentum
         this->dirichlet_function_strings[2] = "if(x<0.5,2.5,0.25)";  // BC for energy
         this->use_exact_solution_as_BC = false;
         // initial conditions
         this->initial_conditions_strings[0] = "if(x<0.5,1.0,0.125)"; // IC for density
         this->initial_conditions_strings[1] = "0";                   // IC for x-momentum
         this->initial_conditions_strings[2] = "if(x<0.5,2.5,0.25)";  // IC for energy
         // exact solution
         this->has_exact_solution = false;
         // physical constants
         gamma = 1.4;
         break;
      }
      default:
         Assert(false,ExcNotImplemented());
   }
}

/** \fn void Euler<dim>::compute_cell_ss_residual()
 *  \brief Computes the contribution of the steady-state residual
 *         from the current cell.
 */
template <int dim>
void Euler<dim>::compute_cell_ss_residual(FEValues<dim> &fe_values,
                                          const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          Vector<double> &cell_residual)
{
   // reinitialize fe values for cell
   fe_values.reinit(cell);

   // get current solution values and gradients
   std::vector<double>         density          (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > density_gradient (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > momentum         (this->n_q_points_cell);
   //std::vector<Tensor<dim,dim> > momentum_symmetric_gradient(this->n_q_points_cell);
   std::vector<double>         energy  (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > velocity(this->n_q_points_cell);
   std::vector<double>         pressure(this->n_q_points_cell);
   fe_values[density_extractor].get_function_values    (this->current_solution, density);
   fe_values[density_extractor].get_function_gradients (this->current_solution, density_gradient);
   fe_values[momentum_extractor].get_function_values   (this->current_solution, momentum);
   //fe_values[momentum_extractor].get_function_symmetric_gradients (this->current_solution, momentum_symmetric_gradient);
   fe_values[energy_extractor].get_function_values     (this->current_solution, energy);

   // compute velocity and pressure
   // loop over quadrature points
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
   {
      Assert(density[q] > 1e-30,ExcDivideByZero());
      velocity[q] = momentum[q] / density[q];
      pressure[q] = compute_pressure(momentum[q],energy[q]);
   }

   // create identity tensor
   SymmetricTensor<2,dim> identity_tensor = unit_symmetric_tensor<dim>();

   // allocate memory for intermediate tensors
   Tensor<2,dim> velocity_times_momentum;
   Tensor<2,dim> velocity_times_momentum_plus_pressure;

   // loop over test functions
   for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      // loop over quadrature points
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      {
         // compute density term
         double density_term = fe_values[density_extractor].gradient(i,q) * momentum[q];

         // compute momentum term
         outer_product(velocity_times_momentum, velocity[q], momentum[q]);
         velocity_times_momentum_plus_pressure = velocity_times_momentum + pressure[q]*identity_tensor;
         double momentum_term = double_contract(velocity_times_momentum_plus_pressure, fe_values[momentum_extractor].gradient(i,q));

         // compute energy term
         double energy_term = fe_values[energy_extractor].gradient(i,q) * velocity[q]*(energy[q] + pressure[q]);

         // sum up terms into cell residual
         cell_residual(i) -= ((
                                density_term
                                  //-density_viscosity[q]*density_gradient[q]
                             +  
                                momentum_term
                                  //-viscosity[q]*velocity_symmetric_gradient[q]
                             +
                                energy_term

                             ) * fe_values.JxW(q));
      }
}

/** \fn void Euler<dim>::compute_face_ss_residual()
 *  \brief Computes the contribution of the steady-state residual
 *         from the faces of the current cell.
 */
template <int dim>
void Euler<dim>::compute_face_ss_residual(FEFaceValues<dim> &fe_face_values,
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
         switch (euler_parameters.viscosity_type)
         {
            case EulerParameters<dim>::constant:
            {
               for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                  viscosity_face(q) = euler_parameters.constant_viscosity_value;
               break;
            }
            case EulerParameters<dim>::first_order:
            {
               // get max velocity on cell
               std::vector<double> local_solution(this->n_q_points_face);
               fe_face_values.get_function_values(this->current_solution, local_solution);
               double max_velocity = 0.0;
               for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                  max_velocity = std::max( max_velocity, local_solution[q]);

               // compute first-order viscosity
               double cell_diameter = cell->diameter();
               double viscosity_value = 0.5 * euler_parameters.first_order_viscosity_coef * cell_diameter * max_velocity;
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
*/
}

/** \fn double Euler<dim>::compute_pressure(const Tensor<1,dim> &momentum,
 *                                          const double        &energy) const
 *  \brief Computes the pressure at a quadrature point.
 */
template <int dim>
double Euler<dim>::compute_pressure(const Tensor<1,dim> &momentum,
                                    const double        &energy) const
{
   return (gamma - 1.0)*(energy - 0.5*std::pow(momentum.norm(),2));
}

/** \fn Tensor<1,dim> Euler<dim>::flux_derivative(const double u)
 *  \brief Computes the derivative of the flux function with respect to
 *         the solution vector.
 *  \param u solution at a point
 *  \return derivative of the flux with respect to u
 */
template <int dim>
Tensor<1,dim> Euler<dim>::flux_derivative(const double u)
{
   Tensor<1,dim> dfdu;
   Assert(false,ExcNotImplemented());
   
/*
   for (unsigned int d = 0; d < dim; ++d)
      dfdu[d] = u;
*/
   return dfdu;
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
void Euler<dim>::compute_entropy(const Vector<double> &solution,
                                 FEValues<dim>        &fe_values,
                                 Vector<double>       &entropy) const
{
Assert(false,ExcNotImplemented());
/*
   std::vector<double> velocity(this->n_q_points_cell);
   fe_values.get_function_values(solution, velocity);

   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      entropy(q) = 0.5*std::pow(velocity[q],2);

   return density/(euler_parameters.gamma)
          * std::log(pressure/std::pow(density,euler_parameters.gamma));
*/
}

/** \fn void Burgers<dim>::compute_entropy_derivative(const Vector<double> &solution,
 *                                                    FEValues<dim> &fe_values,
 *                                                    Vector<double> &entropy_derivative) const
 *  \brief Computes entropy derivative at each quadrature point in cell
 *  \param solution solution
 *  \param fe_values FEValues object
 *  \param entropy_derivative entropy derivative values at each quadrature point in cell
 */
template <int dim>
void Euler<dim>::compute_entropy_derivative(const Vector<double> &solution,
                                            FEValues<dim>        &fe_values,
                                            Vector<double>       &entropy_derivative) const
{
Assert(false,ExcNotImplemented());
/*
   std::vector<double> velocity(this->n_q_points_cell);
   fe_values[momentum].get_function_values(solution, momentum);

   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      entropy_derivative(q) = velocity[q];
*/
}

/** \fn void Euler<dim>::output_solution() const
 *  \brief Outputs the solution to .gpl if 1-D and otherwise to .vtk.
 */
template <int dim>
void Euler<dim>::output_solution () const
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
   
      // number used in output file name for different times
      static unsigned int output_file_number = 0;

      if (dim == 1)
      {
         std::string filename = "output/solution-" +
                                Utilities::int_to_string (output_file_number, 3) +
                                ".gpl";
         std::ofstream output (filename.c_str());

         unsigned int n_cells = this->triangulation.n_active_cells();
         unsigned int total_n_q_points = n_cells * this->n_q_points_cell;
         
         // profiles for each solution variable
         std::vector<std::pair<double,double> > density_profile  (total_n_q_points);
         std::vector<std::pair<double,double> > momentum_profile (total_n_q_points);
         std::vector<std::pair<double,double> > energy_profile   (total_n_q_points);
   
         // solution variables at each quadrature point in cell
         std::vector<double>         density  (this->n_q_points_cell);
         std::vector<Tensor<1,dim> > momentum (this->n_q_points_cell);
         std::vector<double>         energy   (this->n_q_points_cell);
   
         // loop over cells
         FEValues<dim> fe_values (this->fe, this->cell_quadrature, update_values | update_quadrature_points);
         unsigned int i = 0;
         typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
                                                        endc = this->dof_handler.end();
         for (; cell != endc; ++cell)
         {
            fe_values.reinit(cell);
            fe_values[density_extractor].get_function_values  (this->current_solution, density);
            fe_values[momentum_extractor].get_function_values (this->current_solution, momentum);
            fe_values[energy_extractor].get_function_values   (this->current_solution, energy);
   
            for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
            {
               // pair the solution with its x-value
               const Point<dim> q_point = fe_values.quadrature_point(q);
               density_profile[i]  = std::make_pair(q_point(0), density[q]);
               momentum_profile[i] = std::make_pair(q_point(0), momentum[q][0]);
               energy_profile[i]   = std::make_pair(q_point(0), energy[q]);
               ++i;
            }
         }
   
         // sort data by quadrature point
         std::sort(density_profile.begin(),  density_profile.end());
         std::sort(momentum_profile.begin(), momentum_profile.end());
         std::sort(energy_profile.begin(),   energy_profile.end());
   
         // output vector to file
         std::ofstream density_output;
         std::ofstream momentum_output;
         std::ofstream energy_output;
         std::string density_file  = "output/density-"  + Utilities::int_to_string (output_file_number, 3) + ".csv";
         std::string momentum_file = "output/momentum-" + Utilities::int_to_string (output_file_number, 3) + ".csv";
         std::string energy_file   = "output/energy-"   + Utilities::int_to_string (output_file_number, 3) + ".csv";
         density_output.open  (density_file.c_str(),  std::ios::out);
         momentum_output.open (momentum_file.c_str(), std::ios::out);
         energy_output.open   (energy_file.c_str(),   std::ios::out);
         for (i = 0; i < total_n_q_points; ++i)
         {
            density_output  << density_profile[i].first  << "," << density_profile[i].second  << std::endl;
            momentum_output << momentum_profile[i].first << "," << momentum_profile[i].second << std::endl;
            energy_output   << energy_profile[i].first   << "," << energy_profile[i].second   << std::endl;
         }
         density_output.close();
         momentum_output.close();
         energy_output.close();

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
