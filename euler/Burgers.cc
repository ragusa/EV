/**
 * \file Burgers.cc
 * \brief Provides the function definitions for the Burgers class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] params Burgers equation parameters
 */
template <int dim>
Burgers<dim>::Burgers(const BurgersParameters<dim> &params):
   ConservationLaw<dim>(params),
   burgers_parameters(params),
   velocity_extractor(0)
{
} 

template <int dim>
std::vector<std::string> Burgers<dim>::get_component_names ()
{
   std::vector<std::string> names(1,"velocity");
   return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
   Burgers<dim>::get_component_interpretations ()
{
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (1, DataComponentInterpretation::component_is_scalar);

   return data_component_interpretation;
} 

template <int dim>
void Burgers<dim>::define_problem()
{
   switch (burgers_parameters.problem_id)
   {
      case 0: // 1-D, Dirichlet boundary conditions, sin(2*pi*x)
      {
         Assert(dim==1,ExcImpossibleInDim(dim));

         // name of problem
         this->problem_name = "burgers_sin";

         // domain
         double domain_start = 0;
         double domain_width = 1.0;
         this->domain_volume = std::pow(domain_width,dim);
         GridGenerator::hyper_cube(this->triangulation, domain_start,
           domain_start + domain_width);

         // only 1 type of BC: zero Dirichlet; leave boundary indicators as zero
         this->n_boundaries = 1;
         typename Triangulation<dim>::cell_iterator cell = this->triangulation.begin(),
                                                    endc = this->triangulation.end();
         for (; cell != endc; ++cell)
            for (unsigned int face = 0; face < this->faces_per_cell; ++face)
               if (cell->face(face)->at_boundary())
                  cell->face(face)->set_boundary_indicator(0);
         this->boundary_types.resize(this->n_boundaries);
         this->boundary_types[0].resize(this->n_components);
         this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet;
         this->dirichlet_function_strings.resize(this->n_boundaries);
         this->dirichlet_function_strings[0].resize(this->n_components);
         this->dirichlet_function_strings[0][0] = "0";
         this->use_exact_solution_as_BC = false;

         // initial conditions
         this->initial_conditions_strings[0] = "sin(2*pi*x)";

         // exact solution
         this->has_exact_solution = false;

         break;
      }
      case 1: // 1-D Riemann problem, shock wave (u_left > u_right)
      {
         Assert(dim==1,ExcImpossibleInDim(dim));

         // name of problem
         this->problem_name = "burgers_shock";

         // domain
         double domain_start = -1.0;
         double domain_width = 2.0;
         this->domain_volume = std::pow(domain_width,dim);
         GridGenerator::hyper_cube(this->triangulation, domain_start,
           domain_start + domain_width);

         // only 1 type of BC: zero Dirichlet; leave boundary indicators as zero
         this->n_boundaries = 1;
         typename Triangulation<dim>::cell_iterator cell = this->triangulation.begin(),
                                                    endc = this->triangulation.end();
         for (; cell != endc; ++cell)
            for (unsigned int face = 0; face < this->faces_per_cell; ++face)
               if (cell->face(face)->at_boundary())
                  cell->face(face)->set_boundary_indicator(0);
         this->boundary_types.resize(this->n_boundaries);
         this->boundary_types[0].resize(this->n_components);
         this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet;
         this->dirichlet_function_strings.resize(this->n_boundaries);
         this->dirichlet_function_strings[0].resize(this->n_components);
         this->dirichlet_function_strings[0][0] = "if(x<0,1,0)";
         this->use_exact_solution_as_BC = false;

         // initial conditions
         this->initial_conditions_strings[0] = "if(x<0,1,0)";

         // exact solution
         this->has_exact_solution = true;
         this->exact_solution_strings[0] =  "if(x-0.5*t<0,1,0)";

         // create and initialize function parser for exact solution
         std::shared_ptr<FunctionParser<dim> > exact_solution_function_derived =
           std::make_shared<FunctionParser<dim> >(this->parameters.n_components);
         std::map<std::string,double> constants;
         exact_solution_function_derived->initialize("x,t",
                                                     this->exact_solution_strings,
                                                     constants,
                                                     true);
         // point base class pointer to derived class function object
         this->exact_solution_function = exact_solution_function_derived;

         break;
      }
      case 2: // 1-D Riemann problem, rarefaction wave (u_left < u_right)
      {
         Assert(dim==1,ExcImpossibleInDim(dim));

         // name of problem
         this->problem_name = "burgers_rarefaction";

         // domain
         double domain_start = -1.0;
         double domain_width = 2.0;
         this->domain_volume = std::pow(domain_width,dim);
         GridGenerator::hyper_cube(this->triangulation, domain_start,
           domain_start + domain_width);

         // only 1 type of BC: zero Dirichlet; leave boundary indicators as zero
         this->n_boundaries = 1;
         typename Triangulation<dim>::cell_iterator cell = this->triangulation.begin(),
                                                    endc = this->triangulation.end();
         for (; cell != endc; ++cell)
            for (unsigned int face = 0; face < this->faces_per_cell; ++face)
               if (cell->face(face)->at_boundary())
                  cell->face(face)->set_boundary_indicator(0);
         this->boundary_types.resize(this->n_boundaries);
         this->boundary_types[0].resize(this->n_components);
         this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet;
         this->dirichlet_function_strings.resize(this->n_boundaries);
         this->dirichlet_function_strings[0].resize(this->n_components);
         this->dirichlet_function_strings[0][0] = "if(x<0,0,1)";
         this->use_exact_solution_as_BC = false;

         // initial conditions
         this->initial_conditions_strings[0] = "if(x<0,0,1)";

         // exact solution
         this->has_exact_solution = true;
         this->exact_solution_strings[0] = "if(t>0,if(x/t<0,0,if(x/t<1,x/t,1)),if(x<0,0,1))";

         // create and initialize function parser for exact solution
         std::shared_ptr<FunctionParser<dim> > exact_solution_function_derived =
           std::make_shared<FunctionParser<dim> >(this->parameters.n_components);
         std::map<std::string,double> constants;
         exact_solution_function_derived->initialize("x,t",
                                                     this->exact_solution_strings,
                                                     constants,
                                                     true);
         // point base class pointer to derived class function object
         this->exact_solution_function = exact_solution_function_derived;

         break;
      }
      case 3: // Guermond 2-d test problem
      {
         Assert(dim==2,ExcImpossibleInDim(dim));

         // name of problem
         this->problem_name = "burgers_2d";

         // domain
         this->domain_volume = 1.0;
         GridIn<dim> input_grid;
         input_grid.attach_triangulation(this->triangulation);
         std::ifstream input_file("mesh/unit_square.msh");
         input_grid.read_msh(input_file);

         // only 1 type of BC: Dirichlet with exact solution
         this->n_boundaries = 1;
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

         // create and initialize function parser for exact solution
         std::shared_ptr<FunctionParser<dim> > exact_solution_function_derived =
           std::make_shared<FunctionParser<dim> >(this->parameters.n_components);
         std::map<std::string,double> constants;
         exact_solution_function_derived->initialize("x,y,t",
                                                     this->exact_solution_strings,
                                                     constants,
                                                     true);
         // point base class pointer to derived class function object
         this->exact_solution_function = exact_solution_function_derived;

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
void Burgers<dim>::assemble_lumped_mass_matrix()
{
   FEValues<dim> fe_values (this->fe, this->cell_quadrature,
     update_values | update_JxW_values);

   std::vector<types::global_dof_index> local_dof_indices (this->dofs_per_cell);
   FullMatrix<double> local_mass (this->dofs_per_cell, this->dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
                                                  endc = this->dof_handler.end();
   for (; cell != endc; ++cell)
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      local_mass = 0.0;

      // compute local contribution
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
         for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
            for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
            {
               local_mass(i,i) += fe_values[velocity_extractor].value(i,q)
                                  *fe_values[velocity_extractor].value(j,q)
                                  *fe_values.JxW(q);
            }

      // add to global mass matrix with contraints
      this->constraints.distribute_local_to_global (local_mass,
        local_dof_indices, this->lumped_mass_matrix);
   }
}

/**
 * \brief Computes the steady-state residual.
 *
 * Rearranging the Burgers equation,
 * \f[
 *   \frac{\partial u}{\partial t}
 *   = - u\mathbf{v}\cdot\nabla u .
 * \f]
 * Substituting the approximate FEM solution and testing with a test function
 * \f$\varphi_i\f$ gives the weak form for degree of freedom \f$i\f$:
 * \f[
 *   \left(\varphi_i,\frac{\partial u_h}{\partial t}\right)_\Omega
 *   = - \left(\varphi_i,u_h\mathbf{v}\cdot\nabla u_h\right)_\Omega .
 * \f]
 * Adding a viscous bilinear form,
 * \f[
 *   \left(\varphi_i,\frac{\partial u_h}{\partial t}\right)_\Omega
 *   = - \left(\varphi_i,u_h\mathbf{v}\cdot\nabla u_h\right)_\Omega
 *   - \sum\limits_{K\subset S_i}\nu_K\sum\limits_j
 *   U_j b_K(\varphi_i, \varphi_j) .
 * \f]
 * This yields a discrete system
 * \f[
 *   \mathbf{M}\frac{d\mathbf{U}}{dt} = \mathbf{r} ,
 * \f]
 * where \f$\mathbf{M}\f$ is the mass matrix and the steady-state residual
 * \f$\mathbf{r}\f$ is given by
 * \f[
 *   r_i = - \left(\varphi_i,u_h\mathbf{v}\cdot\nabla u_h\right)_\Omega
 *   - \sum\limits_{K\subset S_i}\nu_K\sum\limits_j
 *   U_j b_K(\varphi_i, \varphi_j) .
 * \f]
 *
 *  \param[out] r steady-state residual \f$\mathbf{r}\f$
 */
template <int dim>
void Burgers<dim>::compute_ss_residual(Vector<double> &f)
{
   // reset vector
   f = 0.0;

   FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values |
     update_gradients | update_JxW_values);

   Vector<double> cell_residual(this->dofs_per_cell);
   std::vector<unsigned int> local_dof_indices (this->dofs_per_cell);
   std::vector<double>          solution_values   (this->n_q_points_cell);
   std::vector<Tensor<1, dim> > solution_gradients(this->n_q_points_cell);
   std::vector<Tensor<1, dim> > dfdu              (this->n_q_points_cell);

   //============================================================================
   // inviscid terms
   //============================================================================
   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator
     cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      // reset cell residual
      cell_residual = 0;

      // reinitialize fe values for cell
      fe_values.reinit(cell);
   
      // get current solution values and gradients
      fe_values[velocity_extractor].get_function_values(
        this->new_solution,solution_values);
      fe_values[velocity_extractor].get_function_gradients(
        this->new_solution,solution_gradients);
      
      // loop over quadrature points
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q) {
         // compute derivative of flux
         for (int d = 0; d < dim; ++d)
            dfdu[q][d] = solution_values[q];
         // loop over test functions
         for (unsigned int i = 0; i < this->dofs_per_cell; ++i) {
            cell_residual(i) +=   -fe_values[velocity_extractor].value(i,q)
                                   *dfdu[q]
                                   *solution_gradients[q]
                                   * fe_values.JxW(q);
         }
      }

      // aggregate local residual into global residual
      cell->get_dof_indices(local_dof_indices);
      this->constraints.distribute_local_to_global(cell_residual,
        local_dof_indices, f);
   } // end cell loop

   //============================================================================
   // viscous terms
   //============================================================================
   // if using maximum-principle preserving artificial viscosity, add its
   // bilinear form else use the usual viscous flux contribution
   if (this->parameters.viscosity_type == ConservationLawParameters<dim>::max_principle) {
      this->add_maximum_principle_viscosity_bilinear_form(f);
   } else {
      // loop over cells
      for (cell = this->dof_handler.begin_active(); cell!=endc; ++cell)
      {
         // reset cell residual
         cell_residual = 0;
   
         // reinitialize fe values for cell
         fe_values.reinit(cell);
      
         // get current solution gradients
         fe_values[velocity_extractor].get_function_gradients(this->new_solution,
           solution_gradients);
         
         // loop over quadrature points
         for (unsigned int q = 0; q < this->n_q_points_cell; ++q) {
            // loop over test functions
            for (unsigned int i = 0; i < this->dofs_per_cell; ++i) {
               cell_residual(i) +=   -fe_values[velocity_extractor].gradient(i,q)
                                      *this->viscosity_cell_q[cell](q)
                                      *solution_gradients[q]
                                      * fe_values.JxW(q);
            }
         }

         // aggregate local residual into global residual
         cell->get_dof_indices(local_dof_indices);
         this->constraints.distribute_local_to_global(cell_residual, local_dof_indices, f);
      } // end cell loop
   }
}

/** \brief Computes the contribution of the steady-state residual
 *         from the faces of the current cell.
 *  \param fe_face_values FEFaceValues object.
 *  \param cell current cell.
 *  \param cell_residual residual contribution for the current cell, to be aggregated into the global residual.
 */
/*
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
         fe_face_values[velocity_extractor].get_function_gradients   (this->new_solution,solution_gradients_face);

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
               fe_face_values.get_function_values(this->new_solution, local_solution);
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
*/

/** \brief Computes the steady-state Jacobian matrix and stores in system_matrix.
 */
/*
template <int dim>
void Burgers<dim>::compute_ss_jacobian()
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
      fe_values[velocity_extractor].get_function_values (this->new_solution, velocity);

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
      fe_values[velocity_extractor].get_function_values(this->new_solution, velocity);

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

/** \brief Outputs the solution.
 *  \param[in] time the current time; used in exact solution function
 */
template <int dim>
void Burgers<dim>::output_solution(double time)
{
   DataOut<dim> data_out;
   data_out.attach_dof_handler(this->dof_handler);
   data_out.add_data_vector(this->new_solution,
                            this->component_names,
                            DataOut<dim>::type_dof_data,
                            this->component_interpretations);
   data_out.build_patches ();

   static unsigned int output_file_number = 0;

   if (dim == 1)
   {
      std::string filename = "output/burgers-" +
                             Utilities::int_to_string (output_file_number, 3) +
                             ".gpl";
      std::ofstream output(filename.c_str());
      data_out.write_gnuplot(output);
   }
   else
   {
      std::string filename = "output/burgers-" +
                             Utilities::int_to_string (output_file_number, 3) +
                             ".vtk";
      std::ofstream output(filename.c_str());
      data_out.write_vtk(output);
   }

   // output exact solution if requested
   if (this->parameters.output_exact_solution and
      this->has_exact_solution)
   {
      // set time in exact solution function
      this->exact_solution_function->set_time(time);
      // compute exact solution
      VectorTools::interpolate(
        this->dof_handler,
        *(this->exact_solution_function),
        this->exact_solution);

      DataOut<dim> data_out_exact;
      data_out_exact.attach_dof_handler (this->dof_handler);
      data_out_exact.add_data_vector (this->exact_solution,
                                      this->component_names,
                                      DataOut<dim>::type_dof_data,
                                      this->component_interpretations);
      data_out_exact.build_patches ();

      // if 1-D, output to gnuplot
      if (dim == 1)
      {
         std::string filename_exact = "output/burgers_exact-" +
                                      Utilities::int_to_string (output_file_number, 3) +
                                      ".gpl";
         std::ofstream output_exact(filename_exact.c_str());
         data_out_exact.write_gnuplot(output_exact);
      }
      // else output to vtk
      else
      {
         std::string filename_exact = "output/burgers_exact-" +
                                      Utilities::int_to_string (output_file_number, 3) +
                                      ".vtk";
         std::ofstream output_exact(filename_exact.c_str());
         data_out_exact.write_vtk(output_exact);
      }
   }

   // increment output file number
   ++output_file_number;
}
