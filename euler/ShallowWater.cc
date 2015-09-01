/**
 * \file ShallowWater.cc
 * \brief Provides the function definitions for the ShallowWater class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] params shallow water equation parameters
 */
template <int dim>
ShallowWater<dim>::ShallowWater(const ShallowWaterParameters<dim> &params):
  ConservationLaw<dim>(params),
  burgers_parameters(params),
  height_extractor(0),
  momentum_extractor(1)
{
  // shallow water equations cannot be 3-D
  Assert(dim < 3, ExcImpossibleInDim(dim));
} 

template <int dim>
std::vector<std::string> ShallowWater<dim>::get_component_names ()
{
   std::vector<std::string> names(1 + dim);

   names[0] = "height"
   for (int d = 0; d < dim; ++d) names[1+d] = "momentum";

   return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
   ShallowWater<dim>::get_component_interpretations ()
{
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
     component_interpretations(dim+1);

   component_interpretations[0] =
     DataComponentInterpretation::component_is_scalar;
   for (int d = 0; d < dim; ++d)
     component_interpretations[1+d] =
       DataComponentInterpretation::component_is_part_of_vector;

   return component_interpretations;
} 

template <int dim>
void ShallowWater<dim>::define_problem()
{
   switch (burgers_parameters.problem_id)
   {
      case 0:
      {
         Assert(dim==1,ExcImpossibleInDim(dim));

         // name of problem
         this->problem_name = "dam_break_flat";

         // domain
         double domain_start = -5.0;
         double domain_width = 10.0;
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
         this->boundary_types[0][1] = ConservationLaw<dim>::dirichlet;
         this->dirichlet_function_strings.resize(this->n_boundaries);
         this->dirichlet_function_strings[0].resize(this->n_components);
         this->dirichlet_function_strings[0][0] = "if(x<0,3,1)";
         this->dirichlet_function_strings[0][1] = "0";
         this->use_exact_solution_as_BC = false;

         // initial conditions
         this->initial_conditions_strings[0] = "if(x<0,3,1)";
         this->initial_conditions_strings[1] = "0";

         // exact solution
         this->has_exact_solution = false;
/*
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
*/

         // acceleration due to gravity
         gravity = 1.0;

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
void ShallowWater<dim>::assemble_lumped_mass_matrix()
{
   FEValues<dim> fe_values (this->fe, this->cell_quadrature,
     update_values | update_JxW_values);

   std::vector<types::global_dof_index> local_dof_indices (this->dofs_per_cell);
   FullMatrix<double> local_mass (this->dofs_per_cell, this->dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator
     cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
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
               local_mass(i,i) += (
                 fe_values[height_extractor].value(i,q)
                 * fe_values[height_extractor].value(j,q)
                 + fe_values[momentum_extractor].value(i,q)
                 * fe_values[momentum_extractor].value(j,q)
                 ) * fe_values.JxW(q);
            }

      // add to global mass matrix with contraints
      this->constraints.distribute_local_to_global (local_mass,
        local_dof_indices, this->lumped_mass_matrix);
   }
}

/**
 * \brief Computes the steady-state residual.
 *
 * Rearranging the continuity equation, substituting the approximate FEM
 * solution and testing with a test function \f$\varphi_i^h\f$ gives its
 * weak form for degree of freedom \f$i\f$:
 * \f[
 *   \left(\varphi_i^h,\frac{\partial h_h}{\partial t}\right)_\Omega
 *   = - \left(\varphi_i^h,\nabla\cdot (h \mathbf{u})_h\right)_\Omega .
 * \f]
 * Integrating by parts gives
 * \f[
 *   \left(\varphi_i^h,\frac{\partial h_h}{\partial t}\right)_\Omega
 *   = \left(\nabla\varphi_i^h,(h \mathbf{u})_h\right)_\Omega 
 *   - \left(\varphi_i^h,(h \mathbf{u})_h\cdot\mathbf{n}\right)_{\partial\Omega} .
 * \f]
 * Rearranging the momentum equation, substituting the approximate FEM
 * solution and testing with a test function \f$\varphi_i^{h\mathbf{u}}\f$
 * gives its weak form for degree of freedom \f$i\f$:
 * \f[
 *   \left(\varphi_i^{h\mathbf{u}},\frac{\partial (h\mathbf{u})_h}{\partial t}
 *   \right)_\Omega
 *   = - \left(\varphi_i^{h\mathbf{u}},\nabla\cdot (h\mathbf{u}\otimes\mathbf{u}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\right)_\Omega
 *   + \left(\varphi_i^{h\mathbf{u}},g h_h\nabla b\right)_\Omega .
 * \f]
 * Integrating by parts gives
 * \f[
 *   \left(\varphi_i^{h\mathbf{u}},\frac{\partial (h\mathbf{u})_h}{\partial t}
 *   \right)_\Omega
 *   = \left(\nabla\varphi_i^{h\mathbf{u}},(h\mathbf{u}\otimes\mathbf{u}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\right)_\Omega
 *   - \left(\varphi_i^{h\mathbf{u}},(h\mathbf{u}\otimes\mathbf{u}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\cdot\mathbf{n}\right)_{\partial\Omega}
 *   + \left(\varphi_i^{h\mathbf{u}},g h_h\nabla b\right)_\Omega .
 * \f]
 * This yields a discrete system
 * \f[
 *   \mathbf{M}\frac{d\mathbf{U}}{dt} = \mathbf{r} ,
 * \f]
 * where \f$\mathbf{M}\f$ is the mass matrix and the steady-state residual
 * \f$\mathbf{r}\f$ is given by
 * \f[
 *   r_i = 
 *   \left(\nabla\varphi_i^h,(h \mathbf{u})_h\right)_\Omega 
 *   - \left(\varphi_i^h,(h \mathbf{u})_h\cdot\mathbf{n}\right)_{\partial\Omega}
 *   + \left(\nabla\varphi_i^{h\mathbf{u}},(h\mathbf{u}\otimes\mathbf{u}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\right)_\Omega
 *   - \left(\varphi_i^{h\mathbf{u}},(h\mathbf{u}\otimes\mathbf{u}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\cdot\mathbf{n}\right)_{\partial\Omega}
 *   + \left(\varphi_i^{h\mathbf{u}},g h_h\nabla b\right)_\Omega .
 * \f]
 *
 *  \param[out] r steady-state residual \f$\mathbf{r}\f$
 */
template <int dim>
void ShallowWater<dim>::compute_ss_residual(Vector<double> &f)
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

   //===========================================================================
   // inviscid terms
   //===========================================================================
   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator
     cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      // reset cell residual
      cell_residual = 0;

      // reinitialize fe values for cell
      fe_values.reinit(cell);
   
      // get current solution values
      fe_values[height_extractor].get_function_values(
        this->new_solution, height);
      fe_values[momentum_extractor].get_function_values(
        this->new_solution, momentum);

      // compute inviscid fluxes
      std::vector<Tensor<1,dim> > density_inviscid_flux (this->n_q_points_cell);
      std::vector<Tensor<2,dim> > momentum_inviscid_flux(this->n_q_points_cell);
      compute_inviscid_fluxes(height, momentum,
        height_inviscid_flux, momentum_inviscid_flux);

      // loop over quadrature points
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q) {
         // loop over test functions
         for (unsigned int i = 0; i < this->dofs_per_cell; ++i) {
            cell_residual(i) += (
              // height contributions
              fe_values[height_extractor].gradient(i,q) * height_inviscid_flux[q]
              // momentum contributions
              + double_contract(fe_values[momentum_extractor].gradient(i,q),
              momentum_inviscid_flux[q])
              ) * fe_values.JxW(q);
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
/*
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
*/
}

/**
 * \brief Computes the inviscid fluxes required to be evaluated in cell and face
 *        integrations.
 *
 * \param[in] height  vector of height values
 * \param[in] momentum vector of momentum values
 * \param[out] height_flux  vector of height inviscid flux values
 * \param[out] momentum_flux vector of momentum inviscid flux values
*/
template<int dim>
void ShallowWater<dim>::compute_inviscid_fluxes(
  const std::vector<double>         & height,
  const std::vector<Tensor<1,dim> > & momentum,
  std::vector<Tensor<1,dim> > & height_flux,
  std::vector<Tensor<2,dim> > & momentum_flux
  ) const
{
  // get number of vector elements
  const unsigned int n = density.size();

  // identity tensor
  SymmetricTensor<2,dim> identity_tensor = unit_symmetric_tensor<dim>();

  // compute auxiliary quantities
  std::vector<Tensor<1,dim> > velocity(n);
  compute_velocity(height, momentum, velocity);

  // loop over vector elements
  for (unsigned int q = 0; q < n; ++q)
  {
    // compute density inviscid flux
    height_flux[q] = momentum[q];

    // compute momentum inviscid flux
    Tensor<2,dim> velocity_times_momentum;
    outer_product(velocity_times_momentum, velocity[q], momentum[q]);
    momentum_flux[q] = velocity_times_momentum
      + 0.5*gravity*std::pow(height[q],2)*identity_tensor;
  }
}

/**
 * \brief Computes a vector of velocity values.
 *
 * \param[in] height  vector of height values
 * \param[in] momentum vector of momentum values
 * \param[out] velocity vector of velocity values
 */
template <int dim>
void ShallowWater<dim>::compute_velocity(
   const std::vector<double>         & height,
   const std::vector<Tensor<1,dim> > & momentum,
   std::vector<Tensor<1,dim> >       & velocity) const
{
  unsigned int n = density.size();
  for (unsigned int q = 0; q < n; ++q)
    velocity[q] = momentum[q] / density[q];
}

template <int dim>
void ShallowWater<dim>::update_flux_speeds()
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
void ShallowWater<dim>::compute_entropy(const Vector<double> &solution,
                                   FEValues<dim>        &fe_values,
                                   Vector<double>       &entropy) const
{
   std::vector<double> velocity(this->n_q_points_cell);
   fe_values[velocity_extractor].get_function_values(solution, velocity);

   for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      entropy(q) = 0.5*velocity[q]*velocity[q];
}

template <int dim>
void ShallowWater<dim>::compute_entropy_face(const Vector<double> &solution,
                                        FEFaceValues<dim>    &fe_values_face,
                                        Vector<double>       &entropy) const
{
   std::vector<double> velocity(this->n_q_points_face);
   fe_values_face[velocity_extractor].get_function_values(solution, velocity);

   for (unsigned int q = 0; q < this->n_q_points_face; ++q)
      entropy(q) = 0.5*velocity[q]*velocity[q];
}

template <int dim>
void ShallowWater<dim>::compute_divergence_entropy_flux (const Vector<double> &solution,
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
void ShallowWater<dim>::output_solution(double time)
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
