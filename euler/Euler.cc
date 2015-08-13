/** \file Euler.cc
 *  \brief Provides the function definitions for the Euler class.
 */

/** \brief Constructor for the Euler class.
 *  \param params Euler equation parameters
 */
template <int dim>
Euler<dim>::Euler(const EulerParameters<dim> &params):
  ConservationLaw<dim>(params),
  euler_parameters(params),
  density_extractor(0),
  momentum_extractor(1),
  energy_extractor(1+dim)
{
} 

/** \brief Returns the names of each component.
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

/** \brief Returns the interpretations for each component.
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

/** \brief Create the domain, compute volume, define initial
 *         conditions, and define boundary conditions and exact
 *         solution if it exists.
 */
template <int dim>
void Euler<dim>::define_problem()
{
   switch (euler_parameters.problem_id)
   {
      case 0: // 1-D Sod shock tube problem
      {
         Assert(dim==1,ExcImpossibleInDim(dim));

         // create domain
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

         // set boundary conditions type for each boundary and component
         this->boundary_types.resize(this->n_boundaries);
         this->boundary_types[0].resize(this->n_components);
         this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet; // density has Dirichlet BC
         this->boundary_types[0][1] = ConservationLaw<dim>::dirichlet; // x-momentum has Dirichlet BC
         this->boundary_types[0][2] = ConservationLaw<dim>::dirichlet; // energy has Dirichlet BC

         // set function strings to be parsed for dirichlet boundary condition functions
         this->use_exact_solution_as_BC = false;
         this->dirichlet_function_strings.resize(this->n_boundaries);
         for (unsigned int boundary = 0; boundary < this->n_boundaries; ++boundary) {
            this->dirichlet_function_strings[boundary].resize(this->n_components);
            this->dirichlet_function_strings[boundary][0] = "if(x<0.5,1.0,0.125)"; // BC for density
            this->dirichlet_function_strings[boundary][1] = "0";                   // BC for x-momentum
            this->dirichlet_function_strings[boundary][2] = "if(x<0.5,2.5,0.25)";  // BC for energy
         }

         // initial conditions for each solution component
         this->initial_conditions_strings[0] = "if(x<0.5,1.0,0.125)";
         this->initial_conditions_strings[1] = "0";
         this->initial_conditions_strings[2] = "if(x<0.5,2.5,0.25)";

         // no exact solution coded, although a Riemann solver could be implemented
         this->has_exact_solution = false;

         // physical constants
         gamma = 1.4;

         break;
      }
      case 1: // 1-D Leblanc tube problem
      {
         Assert(dim==1,ExcImpossibleInDim(dim));

         // create domain
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

         // set boundary conditions type for each boundary and component
         this->boundary_types.resize(this->n_boundaries);
         this->boundary_types[0].resize(this->n_components);
         this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet; // density has Dirichlet BC
         this->boundary_types[0][1] = ConservationLaw<dim>::dirichlet; // x-momentum has Dirichlet BC
         this->boundary_types[0][2] = ConservationLaw<dim>::dirichlet; // energy has Dirichlet BC

         // set function strings to be parsed for dirichlet boundary condition functions
         this->use_exact_solution_as_BC = false;
         this->dirichlet_function_strings.resize(this->n_boundaries);
         for (unsigned int boundary = 0; boundary < this->n_boundaries; ++boundary) {
            this->dirichlet_function_strings[boundary].resize(this->n_components);
            this->dirichlet_function_strings[boundary][0] = "if(x<0.5,1.0,0.001)";   // BC for density
            this->dirichlet_function_strings[boundary][1] = "0";                     // BC for x-momentum
            this->dirichlet_function_strings[boundary][2] = "if(x<0.5,0.1,1.0e-10)"; // BC for energy
         }

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
      case 2: // 2-D Noh Problem
      {
         // this is a 2-D problem
         Assert(dim==2,ExcImpossibleInDim(dim));

         // create domain
         this->domain_volume = 1.0; // domain is the unit hypercube, so domain volume is 1^dim
         GridIn<dim> input_grid;
         input_grid.attach_triangulation(this->triangulation);
         std::ifstream input_file("mesh/unit_square.msh");
         input_grid.read_msh(input_file);

         // four boundaries: each side of unit square
         this->n_boundaries = 4;

         // set boundary indicators
         double small_number = 1.0e-15;
         typename Triangulation<dim>::cell_iterator cell = this->triangulation.begin(),
            endc = this->triangulation.end();
         for (; cell != endc; ++cell)
            for (unsigned int face = 0; face < this->faces_per_cell; ++face)
               if (cell->face(face)->at_boundary()) {
                  Point<dim> face_center = cell->face(face)->center();
                  if (face_center(1) < small_number) {            // y = 0 boundary
                     cell->face(face)->set_boundary_indicator(0);
                  } else if (face_center(0) > 1.0-small_number) { // x = 1 boundary
                     cell->face(face)->set_boundary_indicator(1);
                  } else if (face_center(1) > 1.0-small_number) { // y = 1 boundary
                     cell->face(face)->set_boundary_indicator(2);
                  } else if (face_center(0) < small_number) {     // x = 0 boundary
                     cell->face(face)->set_boundary_indicator(3);
                  } else {
                     // all faces should have satisfied one of the conditions
                     std::cout << "x = " << face_center(0) << std::endl;
                     std::cout << "y = " << face_center(1) << std::endl;
                     Assert(false,ExcInternalError());
                  }
               }
         // set boundary conditions type for each boundary and component
         this->boundary_types.resize(this->n_boundaries);
         this->dirichlet_function_strings.resize(this->n_boundaries);
         for (unsigned int boundary = 0; boundary < this->n_boundaries; ++boundary) {
            this->boundary_types[boundary]            .resize(this->n_components);
            this->dirichlet_function_strings[boundary].resize(this->n_components);
            switch (boundary)
            {
               case 0: // y = 0 boundary
               {
                  // all are reflective (zero Neumann)
                  for (unsigned int component = 0; component < this->n_components; ++component)
                     this->boundary_types[boundary][component] = ConservationLaw<dim>::neumann;
                  break;
               }
               case 1: // x = 1 boundary
               {
                  // all are Dirichlet
                  for (unsigned int component = 0; component < this->n_components; ++component)
                     this->boundary_types[boundary][component] = ConservationLaw<dim>::dirichlet;
                  this->dirichlet_function_strings[boundary][0] = "if(sqrt(x^2+y^2)<t/3.0,16,1)"; // density
                  this->dirichlet_function_strings[boundary][1] = "if(sqrt(x^2+y^2)<t/3.0,0,-x/sqrt(x^2+y^2))"; // mx
                  this->dirichlet_function_strings[boundary][2] = "if(sqrt(x^2+y^2)<t/3.0,0,-y/sqrt(x^2+y^2))"; // my
                  this->dirichlet_function_strings[boundary][3] = "if(sqrt(x^2+y^2)<t/3.0,16.0/3.0/(5.0/3.0-1),";
                  this->dirichlet_function_strings[boundary][3] += "1e-9/(5.0/3.0-1)+0.5)"; // energy
                  break;
               }
               case 2: // y = 1 boundary
               {
                  // all are Dirichlet
                  for (unsigned int component = 0; component < this->n_components; ++component)
                     this->boundary_types[boundary][component] = ConservationLaw<dim>::dirichlet;
                  this->dirichlet_function_strings[boundary][0] = "if(sqrt(x^2+y^2)<t/3.0,16,1)"; // density
                  this->dirichlet_function_strings[boundary][1] = "if(sqrt(x^2+y^2)<t/3.0,0,-x/sqrt(x^2+y^2))"; // mx
                  this->dirichlet_function_strings[boundary][2] = "if(sqrt(x^2+y^2)<t/3.0,0,-y/sqrt(x^2+y^2))"; // my
                  this->dirichlet_function_strings[boundary][3] = "if(sqrt(x^2+y^2)<t/3.0,16.0/3.0/(5.0/3.0-1),";
                  this->dirichlet_function_strings[boundary][3] += "1e-9/(5.0/3.0-1)+0.5)"; // energy
                  break;
               }
               case 3: // x = 0 boundary
               {
                  // all are reflective (zero Neumann)
                  for (unsigned int component = 0; component < this->n_components; ++component)
                     this->boundary_types[boundary][component] = ConservationLaw<dim>::neumann;
                  break;
               }
               default:
               {
                  std::cout << "boundary indicator is " << boundary << std::endl;
                  Assert(false,ExcInternalError());
                  break;
               }
            }
         }
         this->use_exact_solution_as_BC = false;

         // initial conditions for each solution component
         this->initial_conditions_strings[0] = "1";
         this->initial_conditions_strings[1] = "if(x==0,0,-x/sqrt(x^2+y^2))";
         this->initial_conditions_strings[2] = "if(y==0,0,-y/sqrt(x^2+y^2))";
         this->initial_conditions_strings[3] = "1e-9/(5.0/3.0-1)+0.5";

         // exact solution
         this->has_exact_solution = true;
         this->exact_solution_strings[0] = "if(sqrt(x^2+y^2)<t/3.0,16,1)"; // density
         this->exact_solution_strings[1] = "if(sqrt(x^2+y^2)<t/3.0,0,-x/sqrt(x^2+y^2))"; // mx
         this->exact_solution_strings[2] = "if(sqrt(x^2+y^2)<t/3.0,0,-y/sqrt(x^2+y^2))"; // my
         this->exact_solution_strings[3] = "if(sqrt(x^2+y^2)<t/3.0,16.0/3.0/(5.0/3.0-1),";
         this->exact_solution_strings[3] += "1e-9/(5.0/3.0-1)+0.5)"; // energy

         // physical constants
         gamma = 5.0/3.0;

         break;
      }
      default:
      {
         Assert(false,ExcNotImplemented());
         break;
      }
   }
}

template <int dim>
void Euler<dim>::compute_ss_residual(Vector<double> &f)
{
   // reset vector
   f = 0.0;

   // NOTE: may need to add more update flags in this constructor
   FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values
     | update_gradients | update_JxW_values);
   FEFaceValues<dim> fe_face_values(this->fe, this->face_quadrature, update_values
     | update_normal_vectors | update_JxW_values);

   Vector<double> cell_residual(this->dofs_per_cell);
   std::vector<unsigned int> local_dof_indices (this->dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
                                                  endc = this->dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      // reset cell residual
      cell_residual = 0;

      // compute local residual
      compute_cell_ss_residual(fe_values, fe_face_values, cell, cell_residual);

      // aggregate local residual into global residual
      cell->get_dof_indices(local_dof_indices);
      this->constraints.distribute_local_to_global(cell_residual, local_dof_indices, f);
   } // end cell loop
}

/** \brief Computes the contribution of the steady-state residual
 *         from the current cell, including face contributions.
 *  \param[in] fe_values FEValues object
 *  \param[in] fe_face_values FEFaceValues object
 *  \param[in] cell cell iterator for current cell
 *  \param[out] cell_residual cell contribution to global right hand side
 */
template <int dim>
void Euler<dim>::compute_cell_ss_residual(
  FEValues<dim>     & fe_values,
  FEFaceValues<dim> & fe_face_values,
  const typename DoFHandler<dim>::active_cell_iterator & cell,
  Vector<double> & cell_residual)
{
   // reinitialize fe values for cell
   fe_values.reinit(cell);

   // get solution values and gradients
   std::vector<double>         density          (this->n_q_points_cell);
   //std::vector<Tensor<1,dim> > density_gradient (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > momentum         (this->n_q_points_cell);
   //std::vector<Tensor<2,dim> > momentum_gradient(this->n_q_points_cell);
   //std::vector<SymmetricTensor<2,dim> > momentum_symmetric_gradient(
   //  this->n_q_points_cell);
   //std::vector<double>         momentum_divergence (this->n_q_points_cell);
   std::vector<double>         energy           (this->n_q_points_cell);
   //std::vector<Tensor<1,dim> > energy_gradient  (this->n_q_points_cell);
/*
   std::vector<double>         internal_energy  (this->n_q_points_cell);
   std::vector<double>         temperature      (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > velocity         (this->n_q_points_cell);
   std::vector<double>         pressure         (this->n_q_points_cell);
   std::vector<double>         dpdrho           (this->n_q_points_cell);
   std::vector<double>         dpdmx            (this->n_q_points_cell);
   std::vector<double>         dpdE             (this->n_q_points_cell);
*/
   fe_values[density_extractor].get_function_values    (this->new_solution, density);
   //fe_values[density_extractor].get_function_gradients (this->new_solution, density_gradient);
   fe_values[momentum_extractor].get_function_values   (this->new_solution, momentum);
   //fe_values[momentum_extractor].get_function_gradients(this->new_solution, momentum_gradient);
   //fe_values[momentum_extractor].get_function_symmetric_gradients (this->new_solution, momentum_symmetric_gradient);
   //fe_values[momentum_extractor].get_function_divergences (this->new_solution, momentum_divergence);
   fe_values[energy_extractor].get_function_values     (this->new_solution, energy);
   //fe_values[energy_extractor].get_function_gradients  (this->new_solution, energy_gradient);

   // compute velocity and thermodynamic properties
/*
   compute_velocity             (density, momentum, velocity);
   compute_internal_energy_cell (internal_energy, density, momentum, energy);
   compute_temperature_cell     (temperature, internal_energy);
   compute_pressure_cell        (pressure, dpdrho, dpdmx, dpdE, density, momentum, energy);
*/

   // check that density is nonzero
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      Assert(density[q] > 1.0e-30, ExcDivideByZero());

/*
   // create identity tensor
   SymmetricTensor<2,dim> identity_tensor = unit_symmetric_tensor<dim>();

   // allocate memory for intermediate tensors
   Tensor<2,dim> velocity_times_momentum;
   Tensor<2,dim> density_viscous_flux_times_velocity;
   Tensor<2,dim> momentum_times_density_gradient;
*/

   // compute inviscid fluxes
   std::vector<Tensor<1,dim> > density_flux (this->n_q_points_cell);
   std::vector<Tensor<2,dim> > momentum_flux(this->n_q_points_cell);
   std::vector<Tensor<1,dim> > energy_flux  (this->n_q_points_cell);
   compute_inviscid_fluxes(density, momentum, energy,
     density_flux, momentum_flux, energy_flux);

   // loop over quadrature points
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
   {
      // compute symmetric gradient of velocity (only able to query symmetric gradient of momentum)
/*
      Tensor<2,dim> aux2;
      outer_product(aux2,momentum[q],density_gradient[q]);
      Tensor<2,dim> velocity_symmetric_gradient = (momentum_gradient[q]*density[q] -aux2)/(density[q]*density[q]);

      // compute viscous fluxes
      double nu = this->viscosity_cell_q[cell](q);
      Tensor<1,dim> density_viscous_flux = -nu*density_gradient[q];

      outer_product(density_viscous_flux_times_velocity, density_viscous_flux, velocity[q]);
      Tensor<2,dim> momentum_viscous_flux = -nu*density[q]*velocity_symmetric_gradient
         + density_viscous_flux_times_velocity;

      double kappa = nu; //euler_parameters.prandtl / (gamma - 1.0) * mu;
      double sp = -pressure[q]/(temperature[q]*density[q]*density[q]);
      double se = 1.0/temperature[q];
      outer_product(momentum_times_density_gradient, momentum[q], density_gradient[q]);
      double velocity_divergence = (momentum_divergence[q]*density[q] - momentum[q]*density_gradient[q])/(density[q]*density[q]);
      Tensor<1,dim> pressure_gradient = (gamma - 1.0)*(energy_gradient[q]
        - velocity_divergence * momentum[q] - 0.5*velocity[q]*velocity[q]*density_gradient[q]);
      Tensor<1,dim> internal_energy_gradient = (pressure_gradient*density[q]
        - pressure[q]*density_gradient[q]) / ((gamma - 1.0)*density[q]*density[q]);
      Tensor<1,dim> l_flux = (nu - kappa)*density[q]*sp/se*density_gradient[q]
        - nu * internal_energy[q] * density_gradient[q]
        - kappa * density[q] * internal_energy_gradient;
      Tensor<1,dim> h_flux = l_flux - 0.5*velocity[q]*velocity[q]*density_viscous_flux;
      Tensor<1,dim> energy_viscous_flux = h_flux + momentum_viscous_flux * velocity[q];
*/

      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
         // sum up terms into cell residual
         cell_residual(i) += ((
            fe_values[density_extractor].gradient(i,q) *
            //(density_flux + density_viscous_flux)
            density_flux[q]
            + double_contract(fe_values[momentum_extractor].gradient(i,q),
            //(momentum_flux + momentum_viscous_flux))
            momentum_flux[q])
            + fe_values[energy_extractor].gradient(i,q) *
            //(energy_flux + energy_viscous_flux)
            energy_flux[q]
         ) * fe_values.JxW(q));
   }

   // resize flux vectors for face computations
   density_flux .resize(this->n_q_points_face);
   momentum_flux.resize(this->n_q_points_face);
   energy_flux  .resize(this->n_q_points_face);

   // loop over faces
   for (unsigned int face = 0; face < this->faces_per_cell; ++face)
   {
      // add term for boundary faces
      if (cell->at_boundary(face))
      {
        fe_face_values.reinit(cell, face);

        // get values on face
        std::vector<double>         density (this->n_q_points_face);
        std::vector<Tensor<1,dim> > momentum(this->n_q_points_face);
        std::vector<double>         energy  (this->n_q_points_face);
        fe_face_values[density_extractor].get_function_values (this->new_solution, density);
        fe_face_values[momentum_extractor].get_function_values(this->new_solution, momentum);
        fe_face_values[energy_extractor].get_function_values  (this->new_solution, energy);

        // compute inviscid fluxes
        compute_inviscid_fluxes(density, momentum, energy,
          density_flux, momentum_flux, energy_flux);
/*
         std::vector<Tensor<1, dim> > solution_gradients_face(this->n_q_points_face);
         fe_face_values[velocity].get_function_gradients   (this->new_solution,solution_gradients_face);

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
               fe_face_values.get_function_values(this->new_solution, local_solution);
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
*/
         // loop over test functions
         for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
            // loop over quadrature points
            for (unsigned int q = 0; q < this->n_q_points_face; ++q)
               cell_residual(i) += (/*fe_face_values[velocity].value(i,q)
                 * viscosity_face(q)
                 * solution_gradients_face[q]
                 * fe_face_values.normal_vector(q)*/
                 - (fe_face_values[density_extractor].value(i,q)  * density_flux[q]
                  + fe_face_values[momentum_extractor].value(i,q) * momentum_flux[q]
                  + fe_face_values[energy_extractor].value(i,q)   * energy_flux[q]
                 ) * fe_face_values.normal_vector(q)
                 )
                 * fe_face_values.JxW(q);
      }
   }
}

/** \brief Computes the contribution of the steady-state residual
 *         from the faces of the current cell.
 *  \param [in] fe_face_values FEFaceValues object
 *  \param [in] cell cell iterator for current cell
 *  \param [out] cell_residual cell contribution to global right hand side
 */
/*
template <int dim>
void Euler<dim>::compute_face_ss_residual(FEFaceValues<dim> &fe_face_values,
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
         fe_face_values[velocity].get_function_gradients   (this->new_solution,solution_gradients_face);

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
               fe_face_values.get_function_values(this->new_solution, local_solution);
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
}
*/

/**
 * \brief Computes the inviscid fluxes required to be evaluated in cell and face
 *        integrations.
 *
 * \param[in] density  vector of density values
 * \param[in] momentum vector of momentum values
 * \param[in] energy   vector of energy values
 * \param[out] density_flux  vector of density inviscid flux values
 * \param[out] momentum_flux vector of momentum inviscid flux values
 * \param[out] energy_flux   vector of energy inviscid flux values
*/
template<int dim>
void Euler<dim>::compute_inviscid_fluxes(
  const std::vector<double>         & density,
  const std::vector<Tensor<1,dim> > & momentum,
  const std::vector<double>         & energy,
  std::vector<Tensor<1,dim> > & density_flux,
  std::vector<Tensor<2,dim> > & momentum_flux,
  std::vector<Tensor<1,dim> > & energy_flux
  ) const
{
  // get number of vector elements
  const unsigned int n = density.size();

  // identity tensor
  SymmetricTensor<2,dim> identity_tensor = unit_symmetric_tensor<dim>();

  // compute auxiliary quantities
  std::vector<Tensor<1,dim> > velocity(n);
  std::vector<double>         pressure(n);
  compute_velocity(density, momentum, velocity);
  compute_pressure(density, momentum, energy, pressure);

  // loop over vector elements
  for (unsigned int q = 0; q < n; ++q)
  {
    // compute density inviscid flux
    density_flux[q] = momentum[q];

    // compute momentum inviscid flux
    Tensor<2,dim> velocity_times_momentum;
    outer_product(velocity_times_momentum, velocity[q], momentum[q]);
    momentum_flux[q] = velocity_times_momentum + pressure[q]*identity_tensor;

    // compute energy inviscid flux
    energy_flux[q] = velocity[q]*(energy[q] + pressure[q]);
  }
}

/** \brief computes the steady-state Jacobian, which is used if an implicit
 *         time integration scheme is to be used.
 */
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
   std::vector<double>         dpdrho           (this->n_q_points_cell);
   std::vector<double>         dpdmx            (this->n_q_points_cell);
   std::vector<double>         dpdE             (this->n_q_points_cell);

   // derivatives of Euler flux functions
   Tensor<1,dim> dfdu_rho_rho;
   Tensor<1,dim> dfdu_rho_mx;
   Tensor<1,dim> dfdu_rho_E;
   Tensor<1,dim> dfdu_mx_rho;
   Tensor<1,dim> dfdu_mx_mx;
   Tensor<1,dim> dfdu_mx_E;
   Tensor<1,dim> dfdu_E_rho;
   Tensor<1,dim> dfdu_E_mx;
   Tensor<1,dim> dfdu_E_E;

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

      fe_values[density_extractor] .get_function_values   (this->new_solution, density);
      fe_values[density_extractor] .get_function_gradients(this->new_solution, density_gradient);
      fe_values[momentum_extractor].get_function_values   (this->new_solution, momentum);
      fe_values[momentum_extractor].get_function_gradients(this->new_solution, momentum_gradient);
      fe_values[momentum_extractor].get_function_symmetric_gradients (this->new_solution, momentum_symmetric_gradient);
      fe_values[momentum_extractor].get_function_divergences(this->new_solution, momentum_divergence);
      fe_values[energy_extractor]  .get_function_values     (this->new_solution, energy);
      fe_values[energy_extractor]  .get_function_gradients  (this->new_solution, energy_gradient);
      compute_pressure_cell (pressure, dpdrho, dpdmx, dpdE, density, momentum, energy);
      dfdu_rho_rho = 0.0;
      dfdu_rho_mx  = unit_vector_x;
      dfdu_rho_E   = 0.0;
      dfdu_mx_rho  = 0.0;

      // reset cell matrix to zero
      cell_matrix = 0;
      // loop over quadrature points in cell
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      {
         double rho = density[q];
         double mx = momentum[q]*unit_vector_x;
         double E = energy[q];
         double p = pressure[q];

         dfdu_mx_mx[0] = 2.0*mx/rho + dpdmx[q];
         dfdu_mx_E[0]  = dpdE[q];
         dfdu_E_rho[0] = -mx/std::pow(rho,2)*(E+p) + mx/rho*dpdrho[q];
         dfdu_E_mx[0]  = (E+p)/rho + mx/rho*dpdmx[q];
         dfdu_E_E[0]   = dpdE[q];

         for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
         {
            for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
            {
               cell_matrix(i,j) += (
                  // conservation of mass
                  fe_values[density_extractor].gradient(i,q) * dfdu_rho_rho *
                  fe_values[density_extractor].value(j,q)
                  + fe_values[density_extractor].gradient(i,q) * dfdu_rho_mx *
                  //fe_values[momentum_extractor].value(j,q)
                  fe_values.shape_value_component(j,q,1)
                  + fe_values[density_extractor].gradient(i,q) * dfdu_rho_E *
                  fe_values[energy_extractor].value(j,q)

                  // conservation of x-momentum
                  //+ fe_values[momentum_extractor].shape_grad_component(i,q,0) * dfdu_mx_rho *
                  + fe_values.shape_grad_component(i,q,1) * dfdu_mx_rho *
                  fe_values[density_extractor].value(j,q)
                  //+ fe_values[momentum_extractor].shape_grad_component(i,q,0) * dfdu_mx_mx *
                  + fe_values.shape_grad_component(i,q,1) * dfdu_mx_mx *
                  fe_values.shape_value_component(j,q,1)
                  //+ fe_values[momentum_extractor].shape_grad_component(i,q,0) * dfdu_mx_E *
                  + fe_values.shape_grad_component(i,q,1) * dfdu_mx_E *
                  fe_values[energy_extractor].value(j,q)

                  // conservation of energy
                  + fe_values[energy_extractor].gradient(i,q) * dfdu_E_rho *
                  fe_values[density_extractor].value(j,q)
                  + fe_values[energy_extractor].gradient(i,q) * dfdu_E_mx *
                  //fe_values[momentum_extractor].value(j,q)
                  fe_values.shape_value_component(j,q,1)
                  + fe_values[energy_extractor].gradient(i,q) * dfdu_E_E *
                  fe_values[energy_extractor].value(j,q)

               ) * fe_values.JxW(q);
            }
         }
      }
      // get dof indices
      cell->get_dof_indices(local_dof_indices);
      // aggregate cell matrix into global matrix
      this->constraints.distribute_local_to_global(cell_matrix, local_dof_indices, this->system_matrix);
   }
}

/** \brief Computes a vector of velocity values.
 */
template <int dim>
void Euler<dim>::compute_velocity(
   const std::vector<double>         & density,
   const std::vector<Tensor<1,dim> > & momentum,
   std::vector<Tensor<1,dim> >       & velocity) const
{
  unsigned int n = density.size();
  for (unsigned int q = 0; q < n; ++q)
    velocity[q] = momentum[q] / density[q];
}

/** \brief Computes internal energy at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_internal_energy_cell(      std::vector<double>         &internal_energy,
                                              const std::vector<double>         &density,
                                              const std::vector<Tensor<1,dim> > &momentum,
                                              const std::vector<double>         &energy) const
{
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      internal_energy[q] = (energy[q] - 0.5*momentum[q]*momentum[q]/density[q]) / density[q];
}

/** \brief Computes internal energy at each quadrature point in face.
 */
template <int dim>
void Euler<dim>::compute_internal_energy_face(      std::vector<double>         &internal_energy,
                                              const std::vector<double>         &density,
                                              const std::vector<Tensor<1,dim> > &momentum,
                                              const std::vector<double>         &energy) const
{
   for (unsigned int q = 0; q < this->n_q_points_face; ++q)
      internal_energy[q] = (energy[q] - 0.5*momentum[q]*momentum[q]/density[q]) / density[q];
}

/** \brief Computes temperature at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_temperature_cell(      std::vector<double> &temperature,
                                          const std::vector<double> &internal_energy) const
{
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      temperature[q] = (gamma - 1.0)*internal_energy[q];
}

/** \brief Computes temperature at each quadrature point in face.
 */
template <int dim>
void Euler<dim>::compute_temperature_face(      std::vector<double> &temperature,
                                          const std::vector<double> &internal_energy) const
{
   for (unsigned int q = 0; q < this->n_q_points_face; ++q)
      temperature[q] = (gamma - 1.0)*internal_energy[q];
}

/**
 * \brief Computes a vector of pressure values.
 *
 * \param[in] density  density
 * \param[in] momentum momentum
 * \param[in] energy   energy
 * \param[out] pressure pressure
 */
template <int dim>
void Euler<dim>::compute_pressure(
  const std::vector<double>          & density,
  const std::vector<Tensor<1, dim> > & momentum,
  const std::vector<double>          & energy,
  std::vector<double>                & pressure) const
{
   unsigned int n = density.size();

   for (unsigned int q = 0; q < n; ++q)
     pressure[q] = (gamma-1.0) * (energy[q]
       - 0.5*momentum[q]*momentum[q]/density[q]);
}

/** \brief Computes pressure and pressure derivatives at each quadrature point in cell.
 *  \param[out] pressure pressure
 *  \param[out] dpdrho derivative of pressure with respect to density
 *  \param[out] dpdmx derivative of pressure with respect to x-momentum
 *  \param[out] dpdE derivative of pressure with respect to energy
 *  \param[in] density density
 *  \param[in] momentum momentum
 *  \param[in] energy energy
 */
template <int dim>
void Euler<dim>::compute_pressure_cell(      std::vector<double> &pressure,
                                             std::vector<double> &dpdrho,
                                             std::vector<double> &dpdmx,
                                             std::vector<double> &dpdE,
                                       const std::vector<double> &density,
                                       const std::vector<Tensor<1, dim> > &momentum,
                                       const std::vector<double> &energy) const
{
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
   {
      pressure[q] = (gamma-1.0) * (energy[q] - 0.5*momentum[q]*momentum[q]/density[q]);
      dpdrho[q] = 0.0;
      dpdmx[q] = (1.0 - gamma)*momentum[q][0];
      dpdE[q] = gamma - 1.0;
   }
}

/** \brief Computes pressure and pressure derivatives at each quadrature point in face.
 *  \param[out] pressure pressure
 *  \param[out] dpdrho derivative of pressure with respect to density
 *  \param[out] dpdmx derivative of pressure with respect to x-momentum
 *  \param[out] dpdE derivative of pressure with respect to energy
 *  \param[in] density density
 *  \param[in] momentum momentum
 *  \param[in] energy energy
 */
template <int dim>
void Euler<dim>::compute_pressure_face(      std::vector<double> &pressure,
                                             std::vector<double> &dpdrho,
                                             std::vector<double> &dpdmx,
                                             std::vector<double> &dpdE,
                                       const std::vector<double> &density,
                                       const std::vector<Tensor<1, dim> > &momentum,
                                       const std::vector<double> &energy) const
{
   for (unsigned int q = 0; q < this->n_q_points_face; ++q)
   {
      pressure[q] = (gamma-1.0) * (energy[q] - 0.5*momentum[q]*momentum[q]);
      dpdrho[q] = 0.0;
      dpdmx[q] = (1.0 - gamma)*momentum[q][0];
      dpdE[q] = gamma - 1.0;
   }
}

/** \brief Computes speed of sound at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_speed_of_sound(      std::vector<double> &speed_of_sound,
   const std::vector<double> &density,
   const std::vector<double> &pressure) const
   {
   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      speed_of_sound[q] = std::sqrt(gamma*pressure[q]/density[q]);
   }

/** \brief Computes the flux speed at each quadrature point in domain and
 *         finds the max in each cell and the max in the entire domain.
 */
template <int dim>
void Euler<dim>::update_flux_speeds()
{
   FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);
   std::vector<double>         density  (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > momentum (this->n_q_points_cell);
   std::vector<double>         energy   (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > velocity (this->n_q_points_cell);
   std::vector<double>         speed_of_sound (this->n_q_points_cell);
   std::vector<double>         pressure (this->n_q_points_cell);
   std::vector<double>         dpdrho   (this->n_q_points_cell);
   std::vector<double>         dpdmx    (this->n_q_points_cell);
   std::vector<double>         dpdE     (this->n_q_points_cell);
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
      fe_values[density_extractor] .get_function_values(this->new_solution, density);
      fe_values[momentum_extractor].get_function_values(this->new_solution, momentum);
      fe_values[energy_extractor]  .get_function_values(this->new_solution, energy);

      // compute velocity
      compute_velocity(density, momentum, velocity);

      // compute speed of sound
      compute_internal_energy_cell (internal_energy, density, momentum, energy);
      compute_temperature_cell     (temperature, internal_energy);
      compute_pressure_cell        (pressure, dpdrho, dpdmx, dpdE, density, momentum, energy);
      compute_speed_of_sound       (speed_of_sound, density, pressure);

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

/** \brief Computes entropy at each quadrature point in cell
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
   std::vector<Tensor<1,dim> > momentum (this->n_q_points_cell);
   std::vector<double>         energy   (this->n_q_points_cell);
   std::vector<double>         pressure (this->n_q_points_cell);
   std::vector<double>         dpdrho   (this->n_q_points_cell);
   std::vector<double>         dpdmx    (this->n_q_points_cell);
   std::vector<double>         dpdE     (this->n_q_points_cell);
   std::vector<double>         temperature     (this->n_q_points_cell);
   std::vector<double>         internal_energy (this->n_q_points_cell);

   fe_values[density_extractor].get_function_values  (this->new_solution, density);
   fe_values[momentum_extractor].get_function_values (this->new_solution, momentum);
   fe_values[energy_extractor].get_function_values   (this->new_solution, energy);
   compute_internal_energy_cell (internal_energy, density, momentum, energy);
   compute_temperature_cell     (temperature, internal_energy);
   compute_pressure_cell        (pressure, dpdrho, dpdmx, dpdE, density, momentum, energy);

   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      entropy[q] = density[q]/(gamma - 1.0)*std::log(pressure[q]/std::pow(density[q],gamma));
   }

/** \brief Computes entropy at each quadrature point on face
 *  \param solution solution
 *  \param fe_values_face FEFaceValues object
 *  \param entropy entropy values at each quadrature point on face
 */
template <int dim>
void Euler<dim>::compute_entropy_face(const Vector<double> &solution,
   FEFaceValues<dim>    &fe_values_face,
   Vector<double>       &entropy) const
   {
   std::vector<double>         density  (this->n_q_points_face);
   std::vector<Tensor<1,dim> > momentum (this->n_q_points_face);
   std::vector<double>         energy   (this->n_q_points_face);
   std::vector<double> pressure         (this->n_q_points_face);
   std::vector<double> dpdrho           (this->n_q_points_face);
   std::vector<double> dpdmx            (this->n_q_points_face);
   std::vector<double> dpdE             (this->n_q_points_face);
   std::vector<double> temperature      (this->n_q_points_face);
   std::vector<double> internal_energy  (this->n_q_points_face);

   fe_values_face[density_extractor].get_function_values  (this->new_solution, density);
   fe_values_face[momentum_extractor].get_function_values (this->new_solution, momentum);
   fe_values_face[energy_extractor].get_function_values   (this->new_solution, energy);
   compute_internal_energy_face (internal_energy, density, momentum, energy);
   compute_temperature_face     (temperature, internal_energy);
   compute_pressure_face        (pressure, dpdrho, dpdmx, dpdE, density, momentum, energy);

   for (unsigned int q = 0; q < this->n_q_points_face; ++q)
      entropy[q] = density[q]/(gamma - 1.0)*std::log(pressure[q]/std::pow(density[q],gamma));
   }

/**
 * Computes divergence of entropy flux at each quadrature point in cell.
 *
 *  \param solution solution
 *  \param fe_values FEValues object
 *  \param divergence divergence of entropy flux at each quadrature point in cell
 */
template <int dim>
void Euler<dim>::compute_divergence_entropy_flux (const Vector<double> &solution,
                                                  FEValues<dim>        &fe_values,
                                                  Vector<double>       &divergence) const
{
   std::vector<double>         density  (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > momentum (this->n_q_points_cell);
   std::vector<double>         energy   (this->n_q_points_cell);
   std::vector<double>         pressure (this->n_q_points_cell);
   std::vector<double>         dpdrho   (this->n_q_points_cell);
   std::vector<double>         dpdmx    (this->n_q_points_cell);
   std::vector<double>         dpdE     (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > density_gradient    (this->n_q_points_cell);
   std::vector<Tensor<2,dim> > momentum_gradient   (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > energy_gradient     (this->n_q_points_cell);
   std::vector<double>         momentum_divergence (this->n_q_points_cell);
   Tensor<1, dim> pressure_gradient;

   fe_values[density_extractor].get_function_values       (this->new_solution, density);
   fe_values[momentum_extractor].get_function_values      (this->new_solution, momentum);
   fe_values[energy_extractor].get_function_values        (this->new_solution, energy);
   fe_values[density_extractor].get_function_gradients    (this->new_solution, density_gradient);
   fe_values[momentum_extractor].get_function_gradients   (this->new_solution, momentum_gradient);
   fe_values[energy_extractor].get_function_gradients     (this->new_solution, energy_gradient);
   fe_values[momentum_extractor].get_function_divergences (this->new_solution, momentum_divergence);
   compute_pressure_cell (pressure, dpdrho, dpdmx, dpdE, density, momentum, energy);

   for (unsigned int q = 0; q < this->n_q_points_cell; ++q) {
      pressure_gradient = (gamma - 1.0) * energy_gradient[q] - momentum[q] * momentum_gradient[q];
      divergence[q] = momentum_divergence[q] * std::log(pressure[q]/std::pow(density[q],gamma))
         + momentum[q] * (pressure_gradient / pressure[q] - gamma/density[q]*density_gradient[q]);
   }
}

/** \brief Outputs the solution to .gpl if 1-D and otherwise to .vtk.
 *  \param time current time; used in exact solution function
 */
template <int dim>
void Euler<dim>::output_solution (double time)
{
   DataOut<dim> data_out;
   data_out.attach_dof_handler (this->dof_handler);

   data_out.add_data_vector (this->new_solution,
      this->component_names,
      DataOut<dim>::type_dof_data,
      this->component_interpretations);

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
         fe_values[density_extractor].get_function_values  (this->new_solution, density);
         fe_values[momentum_extractor].get_function_values (this->new_solution, momentum);
         fe_values[energy_extractor].get_function_values   (this->new_solution, energy);

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
   
   // output exact solution if requested
   if (this->conservation_law_parameters.output_exact_solution and
      this->has_exact_solution)
   {
      // set time in exact solution function
      this->exact_solution_function.set_time(time);
      // compute exact solution
      VectorTools::interpolate(this->dof_handler,this->exact_solution_function,this->exact_solution);

      DataOut<dim> data_out_exact;
      data_out_exact.attach_dof_handler (this->dof_handler);
      data_out_exact.add_data_vector (this->exact_solution,
                                      this->component_names,
                                      DataOut<dim>::type_dof_data,
                                      this->component_interpretations);
      data_out_exact.build_patches ();

      if (dim == 1)
      {
         Assert(false,ExcNotImplemented());
      }
      else
      {
         std::string filename_exact = "output/euler_exact-" +
                                      Utilities::int_to_string (output_file_number, 3) +
                                      ".vtk";
         std::ofstream output_exact(filename_exact.c_str());
         data_out_exact.write_vtk(output_exact);
      }
   }

   ++output_file_number;
}
