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
      case 1: // 1-D Leblanc tube problem
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
         this->dirichlet_function_strings[0] = "if(x<0.5,1.0,0.001)";   // BC for density
         this->dirichlet_function_strings[1] = "0";                     // BC for x-momentum
         this->dirichlet_function_strings[2] = "if(x<0.5,0.1,1.0e-10)"; // BC for energy
         this->use_exact_solution_as_BC = false;
         // initial conditions
         this->initial_conditions_strings[0] = "if(x<0.5,1.0,0.001)";   // IC for density
         this->initial_conditions_strings[1] = "0";                     // IC for x-momentum
         this->initial_conditions_strings[2] = "if(x<0.5,0.1,1.0e-10)"; // IC for energy
         // exact solution
         this->has_exact_solution = false;
         // physical constants
         gamma = 5.0/3.0;
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
   std::vector<Tensor<2,dim> > momentum_gradient(this->n_q_points_cell);
   std::vector<SymmetricTensor<2,dim> > momentum_symmetric_gradient(this->n_q_points_cell);
   std::vector<double>         momentum_divergence (this->n_q_points_cell);
   std::vector<double>         energy           (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > energy_gradient  (this->n_q_points_cell);
   std::vector<double>         internal_energy  (this->n_q_points_cell);
   std::vector<double>         temperature      (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > velocity         (this->n_q_points_cell);
   std::vector<double>         pressure         (this->n_q_points_cell);
   fe_values[density_extractor].get_function_values    (this->current_solution, density);
   fe_values[density_extractor].get_function_gradients (this->current_solution, density_gradient);
   fe_values[momentum_extractor].get_function_values   (this->current_solution, momentum);
   fe_values[momentum_extractor].get_function_gradients(this->current_solution, momentum_gradient);
   fe_values[momentum_extractor].get_function_symmetric_gradients (this->current_solution, momentum_symmetric_gradient);
   fe_values[momentum_extractor].get_function_divergences (this->current_solution, momentum_divergence);
   fe_values[energy_extractor].get_function_values     (this->current_solution, energy);
   fe_values[energy_extractor].get_function_gradients  (this->current_solution, energy_gradient);

   // compute velocity and thermodynamic properties
   compute_velocity(velocity, density, momentum);
   compute_internal_energy(internal_energy, density, momentum, energy);
   compute_temperature(temperature, internal_energy);
   compute_pressure(pressure, density, temperature);

   // loop over quadrature points
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
   {
      Assert(density[q] > 1e-30,ExcDivideByZero());
   }

   // create identity tensor
   SymmetricTensor<2,dim> identity_tensor = unit_symmetric_tensor<dim>();

   // allocate memory for intermediate tensors
   Tensor<2,dim> velocity_times_momentum;
   Tensor<2,dim> density_viscous_flux_times_velocity;
   Tensor<2,dim> momentum_times_density_gradient;

   // loop over quadrature points
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
   {
      // compute density term
      Tensor<1,dim> density_flux = momentum[q];

      // compute momentum term
      outer_product(velocity_times_momentum, velocity[q], momentum[q]);
      Tensor<2,dim> momentum_flux = velocity_times_momentum + pressure[q]*identity_tensor;

      // compute energy term
      Tensor<1,dim> energy_flux = velocity[q]*(energy[q] + pressure[q]);

      // compute symmetric gradient of velocity (only able to query symmetric gradient of momentum)
      Tensor<2,dim> aux2;
      outer_product(aux2,momentum[q],density_gradient[q]);
      Tensor<2,dim> velocity_symmetric_gradient = (momentum_gradient[q]*density[q] -aux2)/(density[q]*density[q]);

      // compute viscous fluxes
      double nu = this->viscosity_cell_q[cell](q);
      Tensor<1,dim> density_viscous_flux = -nu*density_gradient[q];

      outer_product(density_viscous_flux_times_velocity, density_viscous_flux, velocity[q]);
      Tensor<2,dim> momentum_viscous_flux = -nu*density[q]*velocity_symmetric_gradient
                       + density_viscous_flux_times_velocity;

      double kappa = nu;//euler_parameters.prandtl / (gamma - 1.0) * mu;
      double sp = -pressure[q]/(temperature[q]*density[q]*density[q]);
      double se = 1.0/temperature[q];
      outer_product(momentum_times_density_gradient, momentum[q], density_gradient[q]);
      double velocity_divergence = (momentum_divergence[q]*density[q] - momentum[q]*density_gradient[q])/(density[q]*density[q]);
      Tensor<1,dim> pressure_gradient = (gamma - 1.0)*(energy_gradient[q] - velocity_divergence * momentum[q]
         - 0.5*velocity[q]*velocity[q]*density_gradient[q]);
      Tensor<1,dim> internal_energy_gradient = (pressure_gradient*density[q] - pressure[q]*density_gradient[q])/
         ((gamma - 1.0)*density[q]*density[q]);
      Tensor<1,dim> l_flux = (nu - kappa)*density[q]*sp/se*density_gradient[q]
                       - nu    * internal_energy[q] * density_gradient[q]
                       - kappa * density[q]         * internal_energy_gradient;
      Tensor<1,dim> h_flux = l_flux - 0.5*velocity[q]*velocity[q]*density_viscous_flux;
      Tensor<1,dim> energy_viscous_flux = h_flux + momentum_viscous_flux * velocity[q];

      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
         // sum up terms into cell residual
         cell_residual(i) += ((
                                fe_values[density_extractor].gradient(i,q) *
                                 (density_flux + density_viscous_flux)
                             +  
                                double_contract(fe_values[momentum_extractor].gradient(i,q),
                                 (momentum_flux + momentum_viscous_flux))
                             +
                                fe_values[energy_extractor].gradient(i,q) *
                                 (energy_flux + energy_viscous_flux)
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

template <int dim>
void Euler<dim>::compute_ss_jacobian()
{
	// only been coded in 1-D
	Assert(dim==1,ExcNotImplemented());

	// get current solution values and gradients
	std::vector<double>         density          (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > density_gradient (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > momentum         (this->n_q_points_cell);
	std::vector<Tensor<2,dim> > momentum_gradient(this->n_q_points_cell);
	std::vector<SymmetricTensor<2,dim> > momentum_symmetric_gradient(this->n_q_points_cell);
	std::vector<double>         momentum_divergence (this->n_q_points_cell);
	std::vector<double>         energy           (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > energy_gradient  (this->n_q_points_cell);
	std::vector<double>         internal_energy  (this->n_q_points_cell);
	std::vector<double>         temperature      (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > velocity         (this->n_q_points_cell);
	std::vector<double>         pressure         (this->n_q_points_cell);

	// derivatives of Euler flux functions
	std::vector<Tensor<1,dim> > dfdu_rho_rho (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > dfdu_rho_mx  (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > dfdu_rho_E   (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > dfdu_mx_rho  (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > dfdu_mx_mx   (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > dfdu_mx_E    (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > dfdu_E_rho   (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > dfdu_E_mx    (this->n_q_points_cell);
	std::vector<Tensor<1,dim> > dfdu_E_E     (this->n_q_points_cell);

	// other
	Tensor<1,dim> unit_vector_x; unit_vector_x[0] = 1.0;



	// local DoF indices
	std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);

	// cell matrix
	FullMatrix<double> cell_matrix(this->dofs_per_cell,this->dofs_per_cell);

	// FE values
	FEValues<dim> fe_values(this->fe,this->cell_quadrature,
			update_values | update_gradients | update_JxW_values);

	// ones vector
	Tensor<1,dim> ones_vector;
	for (unsigned int d = 0; d < dim; ++d)
		ones_vector[d] = 1.0;

	// reset steady-state Jacobian to zero
	this->system_matrix = 0.0;
	// loop over cells
	typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
			endc = this->dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		fe_values.reinit(cell);

		fe_values[density_extractor] .get_function_values   (this->current_solution, density);
		fe_values[density_extractor] .get_function_gradients(this->current_solution, density_gradient);
		fe_values[momentum_extractor].get_function_values   (this->current_solution, momentum);
		fe_values[momentum_extractor].get_function_gradients(this->current_solution, momentum_gradient);
		fe_values[momentum_extractor].get_function_symmetric_gradients (this->current_solution, momentum_symmetric_gradient);
		fe_values[momentum_extractor].get_function_divergences(this->current_solution, momentum_divergence);
		fe_values[energy_extractor]  .get_function_values     (this->current_solution, energy);
		fe_values[energy_extractor]  .get_function_gradients  (this->current_solution, energy_gradient);
        dfdu_rho_rho = 0.0;
        dfdu_rho_mx  = unit_vector_x;
        dfdu_rho_E   = 0.0;
        dfdu_mx_rho  = 0.0;//

		// reset cell matrix to zero
		cell_matrix = 0;
		// loop over quadrature points in cell
		for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
			for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
			{
				for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
				{
					cell_matrix(i,j) += (
							fe_values[density_extractor].gradient(i,q) * dfdu_rho_rho(q) *
							fe_values[density_extractor].value(j,q)
							+ fe_values[density_extractor].gradient(i,q) * dfdu_rho_mx(q) *
							fe_values[momentum_extractor].value(j,q)
							+ fe_values[density_extractor].gradient(i,q) * dfdu_rho_E(q) *
							fe_values[energy_extractor].value(j,q)

							+ fe_values[momentum_extractor].gradient(i,q) * dfdu_mx_rho(q) *
							fe_values[density_extractor].value(j,q)
							+ fe_values[momentum_extractor].gradient(i,q) * dfdu_mx_mx(q) *
							fe_values[momentum_extractor].value(j,q)
							+ fe_values[momentum_extractor].gradient(i,q) * dfdu_mx_E(q) *
							fe_values[energy_extractor].value(j,q)

							+ fe_values[energy_extractor].gradient(i,q) * dfdu_E_rho(q) *
							fe_values[density_extractor].value(j,q)
							+ fe_values[energy_extractor].gradient(i,q) * dfdu_E_mx(q) *
							fe_values[momentum_extractor].value(j,q)
							+ fe_values[energy_extractor].gradient(i,q) * dfdu_E_E(q) *
							fe_values[energy_extractor].value(j,q)

							) * fe_values.JxW(q);
				}
			}
		// get dof indices
		cell->get_dof_indices(local_dof_indices);
		// aggregate cell matrix into global matrix
		this->constraints.distribute_local_to_global(cell_matrix, local_dof_indices, this->system_matrix);
	}
}

/** \fn void Euler<dim>::compute_velocity(      std::vector<Tensor<1,dim> > &velocity,
 *                                        const std::vector<double>         &density,
 *                                        const std::vector<Tensor<1,dim> > &momentum) const
 *  \brief Computes velocity at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_velocity(      std::vector<Tensor<1,dim> > &velocity,
                                  const std::vector<double>         &density,
                                  const std::vector<Tensor<1,dim> > &momentum) const
{
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      velocity[q] = momentum[q] / density[q];
}

/** \fn void Euler<dim>::compute_internal_energy(      std::vector<double>         &internal_energy,
 *                                               const std::vector<double>         &density,
 *                                               const std::vector<Tensor<1,dim> > &momentum,
 *                                               const std::vector<double>         &energy) const
 *  \brief Computes internal energy at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_internal_energy(      std::vector<double>         &internal_energy,
                                         const std::vector<double>         &density,
                                         const std::vector<Tensor<1,dim> > &momentum,
                                         const std::vector<double>         &energy) const
{
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      internal_energy[q] = (energy[q] - 0.5*momentum[q]*momentum[q]/density[q]) / density[q];
}

/** \fn void Euler<dim>::compute_temperature(      std::vector<double> &temperature,
 *                                           const std::vector<double> &internal_energy) const
 *  \brief Computes temperature at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_temperature(      std::vector<double> &temperature,
                                     const std::vector<double> &internal_energy) const
{
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      temperature[q] = (gamma - 1.0)*internal_energy[q];
}

/** \fn void Euler<dim>::compute_pressure(      std::vector<double> &pressure,
 *                                        const std::vector<double> &density,
 *                                        const std::vector<double> &temperature) const
 *  \brief Computes pressure at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_pressure(      std::vector<double> &pressure,
                                  const std::vector<double> &density,
                                  const std::vector<double> &temperature) const
{
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      pressure[q] = density[q]*temperature[q];
}

/** \fn void Euler<dim>::compute_speed_of_sound(      std::vector<double> &speed_of_sound,
 *                                              const std::vector<double> &density,
 *                                              const std::vector<double> &pressure) const
 *  \brief Computes speed of sound at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_speed_of_sound(      std::vector<double> &speed_of_sound,
                                        const std::vector<double> &density,
                                        const std::vector<double> &pressure) const
{
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      speed_of_sound[q] = std::sqrt(gamma*pressure[q]/density[q]);
}

/** \fn void Euler<dim>::update_flux_speeds()
 *  \brief Computes the flux speed at each quadrature point in domain and
 *     finds the max in each cell and the max in the entire domain.
 */
template <int dim>
void Euler<dim>::update_flux_speeds()
{
   // This has only been implemented for 1-D: how do I compute eigenvalues for multi-D?
   Assert(dim == 1,ExcNotImplemented());

   FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);
   std::vector<double>         density  (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > momentum (this->n_q_points_cell);
   std::vector<double>         energy   (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > velocity (this->n_q_points_cell);
   std::vector<double>         speed_of_sound (this->n_q_points_cell);
   std::vector<double>         pressure (this->n_q_points_cell);
   std::vector<double>         temperature (this->n_q_points_cell);
   std::vector<double>         internal_energy (this->n_q_points_cell);

   // reset max flux speed
   this->max_flux_speed = 0.0;

   // loop over cells to compute first order viscosity at each quadrature point
   typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
                                                  endc = this->dof_handler.end();
   for (; cell != endc; ++cell)
   {
      fe_values.reinit(cell);
      fe_values[density_extractor].get_function_values (this->current_solution, density);
      fe_values[momentum_extractor].get_function_values(this->current_solution, momentum);
      fe_values[energy_extractor].get_function_values(this->current_solution, energy);
      // compute velocity
      compute_velocity(velocity, density, momentum);
      // compute speed of sound
      compute_internal_energy(internal_energy, density, momentum, energy);
      compute_temperature(temperature, internal_energy);
      compute_pressure(pressure, density, temperature);
      compute_speed_of_sound(speed_of_sound, density, pressure);

      double max_speed_of_sound = 0.0; // max speed of sound on cell
      double max_velocity = 0.0;       // max velocity on cell
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      {
         max_speed_of_sound = std::max( max_speed_of_sound, speed_of_sound[q] );
         max_velocity       = std::max( max_velocity,       velocity[q].norm() );
      }

      // get max flux speed
      this->max_flux_speed_cell[cell] = max_speed_of_sound + max_velocity;
      this->max_flux_speed = std::max(this->max_flux_speed, this->max_flux_speed_cell[cell]);
   }
}

/** \fn void Euler<dim>::compute_entropy(const Vector<double> &solution,
 *                                       FEValues<dim> &fe_values,
 *                                       Vector<double> &entropy) const
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
   std::vector<double>         density  (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > momentum  (this->n_q_points_cell);
   std::vector<double>         energy  (this->n_q_points_cell);
   std::vector<double> pressure  (this->n_q_points_cell);
   std::vector<double> temperature  (this->n_q_points_cell);
   std::vector<double> internal_energy  (this->n_q_points_cell);

   fe_values[density_extractor].get_function_values  (this->current_solution, density);
   fe_values[momentum_extractor].get_function_values (this->current_solution, momentum);
   fe_values[energy_extractor].get_function_values   (this->current_solution, energy);
   compute_internal_energy(internal_energy, density, momentum, energy);
   compute_temperature(temperature, internal_energy);
   compute_pressure(pressure, density, temperature);

   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      entropy[q] = density[q]/(gamma - 1.0)*std::log(pressure[q]/std::pow(density[q],gamma));
}

/** \fn void Euler<dim>::compute_entropy_face(const Vector<double> &solution,
 *                                            FEFaceValues<dim> &fe_values_face,
 *                                            Vector<double> &entropy) const
 *  \brief Computes entropy at each quadrature point on face
 *  \param solution solution
 *  \param fe_values_face FEFaceValues object
 *  \param entropy entropy values at each quadrature point on face
 */
template <int dim>
void Euler<dim>::compute_entropy_face(const Vector<double> &solution,
                                      FEFaceValues<dim>    &fe_values_face,
                                      Vector<double>       &entropy) const
{
Assert(false,ExcNotImplemented());
}

/** \fn void Euler<dim>::compute_divergence_entropy_flux (const Vector<double> &solution,
 *                                                        FEValues<dim> &fe_values,
 *                                                        Vector<double> &divergence) const
 *  \brief Computes divergence of entropy flux at each quadrature point in cell
 *  \param solution solution
 *  \param fe_values FEValues object
 *  \param divergence divergence of entropy flux at each quadrature point in cell
 */
template <int dim>
void Euler<dim>::compute_divergence_entropy_flux (const Vector<double> &solution,
                                                  FEValues<dim>        &fe_values,
                                                  Vector<double>       &divergence) const
{
Assert(false,ExcNotImplemented());
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
