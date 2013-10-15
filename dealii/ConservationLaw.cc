/** \file ConservationLaw.cc
 *  \brief Provides function definitions for the ConservationLaw class.
 */

/** \fn ConservationLaw<dim>::ConservationLaw(const ConservationLawParameters<dim> &params)
 *  \brief Constructor for ConservationLaw class.
 *  \param params conservation law parameters
 */
template <int dim>
ConservationLaw<dim>::ConservationLaw(const ConservationLawParameters<dim> &params):
   conservation_law_parameters(params),
   n_components(params.n_components),
   mapping(),
   fe(FE_Q<dim>(params.degree), params.n_components),
   dofs_per_cell(fe.dofs_per_cell),
   faces_per_cell(GeometryInfo<dim>::faces_per_cell),
   dof_handler(triangulation),
   n_q_points_per_dim(params.n_quadrature_points),
   cell_quadrature(n_q_points_per_dim),
   face_quadrature(n_q_points_per_dim),
   n_q_points_cell(cell_quadrature.size()),
   n_q_points_face(face_quadrature.size()),
   dirichlet_function_strings(params.n_components),
   dirichlet_function        (params.n_components),
   initial_conditions_strings (params.n_components),
   initial_conditions_function(params.n_components),
   exact_solution_strings (params.n_components),
   exact_solution_function(params.n_components)
{}

/** \fn void ConservationLaw<dim>::run()
 *  \brief Runs the entire program.
 *
 *  This function is the uppermost level function that
 *  calls all other functions.
 */
template <int dim>
void ConservationLaw<dim>::run()
{
   // initialize system
   initialize_system();

   // loop over adaptive refinement cycles
   for (unsigned int cycle = 0; cycle < conservation_law_parameters.n_cycle; ++cycle)
   {
      std::cout << std::endl;
      std::cout << "Cycle " << cycle+1 << " of " << conservation_law_parameters.n_cycle << ":" << std::endl;;

      // if in final cycle, set flag to output solution
      if (cycle == conservation_law_parameters.n_cycle-1)
         in_final_cycle = true;
      else
         in_final_cycle = false;

      // adaptively refine mesh if not the first cycle
      if (cycle > 0)
         refine_mesh();

      std::cout << "Number of active cells: ";
      std::cout << triangulation.n_active_cells();
      std::cout << std::endl;

      // setup system; to be applied after each refinement
      setup_system();
   
      // interpolate the initial conditions to the grid
      VectorTools::interpolate(dof_handler,initial_conditions_function,current_solution);
      // apply Dirichlet BC to initial solution or guess
      apply_Dirichlet_BC(0.0);
      // set old solution to the current solution
      old_solution = current_solution;
      // output initial solution
      output_solution();
   
      // begin time stepping loop
      switch (conservation_law_parameters.temporal_integrator)
      {
         case ConservationLawParameters<dim>::runge_kutta: // explicit Runge-Kutta
             solve_runge_kutta();
             break;
         default:
             Assert(false,ExcNotImplemented());
      }
   }

   // output final viscosities if non-constant viscosity used
   switch (conservation_law_parameters.viscosity_type)
   {
      case ConservationLawParameters<dim>::none:
         break;
      case ConservationLawParameters<dim>::constant:
         break;
      case ConservationLawParameters<dim>::first_order:
         output_map(first_order_viscosity_cell_q, "first_order_viscosity");
         break;
      case ConservationLawParameters<dim>::entropy:
         output_map(first_order_viscosity_cell_q, "first_order_viscosity");
         output_map(entropy_viscosity_cell_q, "entropy_viscosity");
         output_map(viscosity_cell_q, "viscosity");
         if (conservation_law_parameters.add_jumps)
            output_map(entropy_viscosity_with_jumps_cell_q, "entropy_viscosity_with_jumps");
         break;
      default:
         Assert(false,ExcNotImplemented());
   }
}

/** \fn void ConservationLaw<dim>::initialize_system()
 *  \brief Initially sets up system. Called only once.
 */
template <int dim>
void ConservationLaw<dim>::initialize_system()
{
   // get component names and interpretations
   component_names           = get_component_names();
   component_interpretations = get_component_interpretations();

   // make grid and refine
   define_problem();
   triangulation.refine_global(conservation_law_parameters.initial_refinement_level);

   // set flag to skip computing face residuals if all BC are Dirichlet
   need_to_compute_face_residual = false;
   bool at_least_one_dirichlet_BC = false;
   for (unsigned int boundary = 0; boundary < n_boundaries; ++boundary)
      for (unsigned int component = 0; component < n_components; ++component)
         if (boundary_types[boundary][component] == dirichlet)
            at_least_one_dirichlet_BC = true;
         else
            need_to_compute_face_residual = true;

   // create constants for parsed functions
   std::map<std::string,double> constants;
   constants["pi"] = numbers::PI;

   // initialize Dirichlet boundary functions if needed
   if (at_least_one_dirichlet_BC && !(use_exact_solution_as_BC))
      dirichlet_function.initialize(FunctionParser<dim>::default_variable_names(),
                                    dirichlet_function_strings,
                                    constants,
                                    true);

   // create exact solution function if there is one
   if (has_exact_solution)
   {
      std::string variables;
      if (dim == 1)
         variables = "x,t";
      else if (dim == 2)
         variables = "x,y,t";
      else if (dim == 3)
         variables = "x,y,z,t";
      else
         Assert(false,ExcInvalidState());

      exact_solution_function.initialize(variables,
                                         exact_solution_strings,
                                         constants,
                                         true);
   }

   // create initial conditions function
   initial_conditions_function.initialize(FunctionParser<dim>::default_variable_names(),
                                          initial_conditions_strings,
                                          constants,
                                          false);

   // initialize Runge-Kutta data if RK is being used
   if (conservation_law_parameters.temporal_integrator ==
      ConservationLawParameters<dim>::runge_kutta)
      initialize_runge_kutta();

}

/** \fn void ConservationLaw<dim>::initialize_runge_kutta()
 *  \brief Assigns Butcher tableau constants.
 */
template <int dim>
void ConservationLaw<dim>::initialize_runge_kutta()
{
   // get RK parameters a, b, and c (Butcher tableau)
   rk.s = 0;
   switch (conservation_law_parameters.runge_kutta_method)
   {
      case ConservationLawParameters<dim>::erk1:
         rk.s = 1;
         break;
      case ConservationLawParameters<dim>::erk2:
         rk.s = 2;
         break;
      case ConservationLawParameters<dim>::erk3:
         rk.s = 3;
         break;
      case ConservationLawParameters<dim>::erk4:
         rk.s = 4;
         break;
      default:
         Assert(false,ExcNotImplemented());
         break;
   }

   // allocate memory for constants
   rk.a.resize(rk.s);
   for (int i = 0; i < rk.s; ++i)
      rk.a[i].reinit(rk.s);
   rk.b.resize(rk.s);
   rk.c.resize(rk.s);

   // assign constants
   switch (conservation_law_parameters.runge_kutta_method)
   {
      case ConservationLawParameters<dim>::erk1:
         rk.b[0] = 1;
         rk.c[0] = 0;
         break;
      case ConservationLawParameters<dim>::erk2:
         rk.a[1][0] = 0.5;
         rk.b[0] = 0;
         rk.b[1] = 1;
         rk.c[0] = 0;
         rk.c[1] = 0.5;
         break;
      case ConservationLawParameters<dim>::erk3:
         rk.a[1][0] = 1.0;
         rk.a[2][0] = 0.25;
         rk.a[2][1] = 0.25;
         rk.b[0] = 1./6;
         rk.b[1] = 1./6;
         rk.b[2] = 4./6;
         rk.c[0] = 0;
         rk.c[1] = 1.0;
         rk.c[2] = 0.5;
         break;
      case ConservationLawParameters<dim>::erk4:
         rk.a[1][0] = 0.5;
         rk.a[2][0] = 0;
         rk.a[2][1] = 0.5;
         rk.a[3][0] = 0;
         rk.a[3][1] = 0;
         rk.a[3][2] = 1;
         rk.b[0] = 1./6;
         rk.b[1] = 1./3;
         rk.b[2] = 1./3;
         rk.b[3] = 1./6;
         rk.c[0] = 0;
         rk.c[1] = 0.5;
         rk.c[2] = 0.5;
         rk.c[3] = 1;
         break;
      default:
         Assert(false,ExcNotImplemented());
         break;
   }

   // allocate vectors for steady-state residual evaluations
   rk.f.resize(rk.s);

   // test to see if last row of A matrix is the same as the b vector;
   // this implies the last stage solution is the new solution and thus
   // the final linear combination is not required.
   rk.solution_computed_in_last_stage = true;
   for (int i = 0; i < rk.s; ++i)
      if (std::abs(rk.a[rk.s-1][i]-rk.b[i]) > 1.0e-30)
         rk.solution_computed_in_last_stage = false;
}

/** \fn void ConservationLaw<dim>::compute_error_for_refinement()
 *  \brief Computes error for adaptive mesh refinement for a time
 *         step and adds it to an error sum for all time steps.
 */
template <int dim>
void ConservationLaw<dim>::compute_error_for_refinement()
{
   Vector<float> estimated_error_per_cell_time_step (triangulation.n_active_cells());

   KellyErrorEstimator<dim>::estimate (dof_handler,
                                       face_quadrature,
                                       // for now, assume no Neumann boundary conditions,
                                       //  so the following argument may be empty
                                       typename FunctionMap<dim>::type(),
                                       current_solution,
                                       estimated_error_per_cell_time_step);

   estimated_error_per_cell += estimated_error_per_cell_time_step;
}

/** \fn void ConservationLaw<dim>::refine_mesh()
 *  \brief Adaptively refines mesh based on a previously computed error per cell.
 */
template <int dim>
void ConservationLaw<dim>::refine_mesh()
{
   GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                    estimated_error_per_cell,
                                                    conservation_law_parameters.refinement_fraction,
                                                    conservation_law_parameters.coarsening_fraction);

   triangulation.execute_coarsening_and_refinement();
}

/** \fn void ConservationLaw<dim>::setup_system()
 *  \brief Sets up the system before solving.
 *
 *  This function is to be applied after each refinement. It
 *  allocates memory, sets up constraints, makes the sparsity pattern,
 *  and reinitializes the system matrix with the sparsity pattern.
 */
template <int dim>
void ConservationLaw<dim>::setup_system ()
{
   // clear maps
   cell_diameter.clear();
   max_flux_speed_cell.clear();
   viscosity_cell_q.clear();
   first_order_viscosity_cell_q.clear();
   entropy_viscosity_cell_q.clear();
   entropy_viscosity_with_jumps_cell_q.clear();
   entropy_cell_q.clear();
   entropy_residual_cell_q.clear();
   max_entropy_residual_cell.clear();
   max_jumps_cell.clear();

   // update cell sizes and minimum cell size
   update_cell_sizes();

   // clear and distribute dofs
   dof_handler.clear();
   dof_handler.distribute_dofs(fe);

   std::cout << "Number of degrees of freedom: ";
   std::cout << dof_handler.n_dofs();
   std::cout << std::endl;

   // make constraints
   constraints.clear();
   DoFTools::make_hanging_node_constraints (dof_handler, constraints);
   for (unsigned int boundary = 0; boundary < n_boundaries; ++boundary)
      for (unsigned int component = 0; component < n_components; ++component)
         if (boundary_types[boundary][component] == dirichlet)
         {
            std::vector<bool> component_mask(n_components, false);
            component_mask[component] = true;
            if (use_exact_solution_as_BC)
            {
               exact_solution_function.set_time(0.0);
               VectorTools::interpolate_boundary_values (dof_handler,
                                                         boundary,
                                                         exact_solution_function,
                                                         constraints,
                                                         component_mask);
            } else {
               dirichlet_function.set_time(0.0);
               VectorTools::interpolate_boundary_values (dof_handler,
                                                         boundary,
                                                         dirichlet_function,
                                                         constraints,
                                                         component_mask);
            }
         }
   constraints.close();

   // create sparsity pattern and compress
   CompressedSparsityPattern compressed_sparsity_pattern (dof_handler.n_dofs() );
   DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern, constraints, false);
   sparsity_pattern.copy_from(compressed_sparsity_pattern);

   // create mass matrix
   mass_matrix.reinit(sparsity_pattern);
   assemble_mass_matrix();

   // resize vectors
   old_solution.reinit(dof_handler.n_dofs());
   current_solution.reinit(dof_handler.n_dofs());
   system_rhs.reinit(dof_handler.n_dofs());
   estimated_error_per_cell.reinit(triangulation.n_active_cells());
   
   // allocate memory for steady-state residual evaluations if Runge-Kutta is used
   if (conservation_law_parameters.temporal_integrator ==
      ConservationLawParameters<dim>::runge_kutta)
      for (int i = 0; i < rk.s; ++i)
         rk.f[i].reinit(dof_handler.n_dofs());

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      viscosity_cell_q[cell]             = Vector<double>(n_q_points_cell);
      first_order_viscosity_cell_q[cell] = Vector<double>(n_q_points_cell);
      entropy_viscosity_cell_q[cell]     = Vector<double>(n_q_points_cell);
      entropy_cell_q[cell]               = Vector<double>(n_q_points_cell);
      entropy_residual_cell_q[cell]      = Vector<double>(n_q_points_cell);
      entropy_viscosity_with_jumps_cell_q[cell]     = Vector<double>(n_q_points_cell);
   }
}

/** \fn void ConservationLaw<dim>::update_cell_sizes()
 *  \brief Updates the cell sizes map and minimum cell size.
 */
template <int dim>
void ConservationLaw<dim>::update_cell_sizes()
{
   // fill cell size map and find minimum cell size
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   // reset minimum cell size to an arbitrary cell size
   minimum_cell_diameter = cell->diameter();

   for (; cell != endc; ++cell)
   {
      cell_diameter[cell] = cell->diameter();
      minimum_cell_diameter = std::min( minimum_cell_diameter, cell_diameter[cell] );
   }
}

/** \fn void ConservationLaw<dim>::assemble_mass_matrix()
 *  \brief Assembles the mass matrix and applies constraints.
 */
template <int dim>
void ConservationLaw<dim>::assemble_mass_matrix ()
{
   mass_matrix = 0.0;

// The below code creates a singular mass matrix for vector-valued problems
/*
   FEValues<dim> fe_values (fe, cell_quadrature, update_values | update_JxW_values);

   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   FullMatrix<double> local_mass (dofs_per_cell, dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      local_mass = 0.0;

      // compute local contribution
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
               local_mass(i,j) +=  fe_values.shape_value(i,q)
                                  *fe_values.shape_value(j,q)
                                  *fe_values.JxW(q);
            }

      // add to global mass matrix with contraints
      constraints.distribute_local_to_global (local_mass, local_dof_indices, mass_matrix);
   }
*/
   Function<dim> *dummy_function = 0;
   MatrixTools::create_mass_matrix( dof_handler,
                                    cell_quadrature,
                                    mass_matrix,
                                    dummy_function,
                                    constraints);

   // output mass matrix
   if (conservation_law_parameters.output_mass_matrix)
   {
      std::ofstream mass_matrix_out ("output/mass_matrix.txt");
      mass_matrix.print_formatted(mass_matrix_out, 10, true, 0, "0", 1);
      mass_matrix_out.close();
   }
}

/** \fn void ConservationLaw<dim>apply_Dirichlet_BC(const double &time)
 *  \brief Applies Dirichlet boundary conditions
 *
 *  This function applies Dirichlet boundary conditions
 *  using the interpolate_boundary_values() tool.
 */
template <int dim>
void ConservationLaw<dim>::apply_Dirichlet_BC(const double &time)
{
   std::map<unsigned int, double> boundary_values;
   for (unsigned int boundary = 0; boundary < n_boundaries; ++boundary)
      for (unsigned int component = 0; component < n_components; ++component)
         if (boundary_types[boundary][component] == dirichlet)
         {
            // create mask to prevent function from being applied to other components
            std::vector<bool> component_mask(n_components, false);
            component_mask[component] = true;
            // fill boundary_values with boundary values
            if (use_exact_solution_as_BC)
            {
               exact_solution_function.set_time(time);
               VectorTools::interpolate_boundary_values (dof_handler,
                                                         boundary,
                                                         exact_solution_function,
                                                         boundary_values,
                                                         component_mask);
            } else {
               dirichlet_function.set_time(time);
               VectorTools::interpolate_boundary_values (dof_handler,
                                                         boundary,
                                                         dirichlet_function,
                                                         boundary_values,
                                                         component_mask);
            }
         }
   // apply boundary values to the solution
   for (std::map<unsigned int, double>::const_iterator it = boundary_values.begin(); it != boundary_values.end(); ++it)
      current_solution(it->first) = (it->second);
}

/** \fn void ConservationLaw<dim>::output_map(std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > &map,
                                              const std::string &output_filename_base)
 *  \brief Outputs a mapped quantity at all quadrature points
 */
template <int dim>
void ConservationLaw<dim>::output_map(std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > &map,
                                      const std::string &output_filename_base)
{
   if (dim == 1)
   {
      // get data from maps into vector of point-data pairs
      unsigned int n_cells = triangulation.n_active_cells();
      unsigned int total_n_q_points = n_cells * n_q_points_cell;
      std::vector<std::pair<double,double> > profile(total_n_q_points);

      FEValues<dim> fe_values (fe, cell_quadrature, update_quadrature_points);
      unsigned int i = 0;
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                     endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
         fe_values.reinit(cell);

         for (unsigned int q = 0; q < n_q_points_cell; ++q)
         {
            const Point<dim> q_point = fe_values.quadrature_point(q);
            profile[i] = std::make_pair(q_point(0), map[cell](q));
            ++i;
         }
      }

      // sort data by quadrature point
      std::sort(profile.begin(), profile.end());

      // output vector to file
      std::ofstream output;
      std::string output_file = "output/" + output_filename_base + ".csv";
      output.open(output_file.c_str(), std::ios::out);
      for (i = 0; i < total_n_q_points; ++i)
         output << profile[i].first << "," << profile[i].second << std::endl;
      output.close();
   }
   else
   {
      std::cout << "Maps were not output because this has only been implemented for 1d" << std::endl;
   }
}

/** \fn void ConservationLaw<dim>::output_map(std::map<typename DoFHandler<dim>::active_cell_iterator, double> &map,
                                              const std::string &output_filename_base)
 *  \brief Outputs a mapped quantity at all quadrature points
 */
template <int dim>
void ConservationLaw<dim>::output_map(std::map<typename DoFHandler<dim>::active_cell_iterator, double> &map,
                                      const std::string &output_filename_base)
{
   if (dim == 1)
   {
      // get data from maps into vector of point-data pairs
      unsigned int n_cells = triangulation.n_active_cells();
      unsigned int total_n_q_points = n_cells * n_q_points_cell;
      std::vector<std::pair<double,double> > profile(total_n_q_points);

      FEValues<dim> fe_values (fe, cell_quadrature, update_quadrature_points);
      unsigned int i = 0;
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                     endc = dof_handler.end();
      for (; cell != endc; ++cell)
      {
         fe_values.reinit(cell);

         for (unsigned int q = 0; q < n_q_points_cell; ++q)
         {
            const Point<dim> q_point = fe_values.quadrature_point(q);
            profile[i] = std::make_pair(q_point(0), map[cell]);
            ++i;
         }
      }

      // sort data by quadrature point
      std::sort(profile.begin(), profile.end());

      // output vector to file
      std::ofstream output;
      std::string output_file = "output/" + output_filename_base + ".csv";
      output.open(output_file.c_str(), std::ios::out);
      for (i = 0; i < total_n_q_points; ++i)
         output << profile[i].first << "," << profile[i].second << std::endl;
      output.close();
   }
   else
   {
      std::cout << "Maps were not output because this has only been implemented for 1d" << std::endl;
   }
}

/** \fn ConservationLaw<dim>::solve_runge_kutta()
 *  \brief Solves the transient using Runge-Kutta.
 *
 *  This function contains the transient loop and solves the
 *  transient using explicit Runge-Kutta:
 *  \f[
 *    \mathbf{M} \mathbf{y}_{n+1} = \mathbf{M} \mathbf{y}_n + h\sum\limits^s_{i=1}b_i \mathbf{f}_i
 *  \f]
 *  where
 *  \f[
 *    \mathbf{f}_i = \mathbf{f}(t_n + c_i h, \mathbf{Y}_i)
 *  \f]
 *  and \f$\mathbf{Y}_i\f$ is computed from the linear solve
 *  \f[
 *    \mathbf{M} \mathbf{Y}_i = \mathbf{M} \mathbf{y}_n + h\sum\limits^{i-1}_{j=1}a_{i,j} \mathbf{f}_i
 *  \f]
 */
template <int dim>
void ConservationLaw<dim>::solve_runge_kutta()
{
   old_time = 0.0;
   unsigned int n = 1; // time step index
   unsigned int next_time_step_output = conservation_law_parameters.output_period;
   double t_end = conservation_law_parameters.final_time;
   bool final_time_not_reached = true;
   while (final_time_not_reached)
   {
      // update max speed for use in CFL computation
      update_flux_speeds();

      // compute dt
      double dt;
      switch (conservation_law_parameters.time_step_size_method)
      {
         case ConservationLawParameters<dim>::constant_dt:
            dt = conservation_law_parameters.time_step_size;
            break;
         case ConservationLawParameters<dim>::cfl_condition:
            dt = compute_dt_from_cfl_condition();
            break;
         default:
            Assert(false,ExcNotImplemented());
            break;
      }
      // check end of transient and shorten last time step if necessary
      if ((old_time+dt) >= t_end)
      {
         dt = t_end - old_time;
         final_time_not_reached = false;
      }
      // compute CFL number
      double cfl = compute_cfl_number(dt);

      std::cout << "   time step " << n << ": t = " << old_time << "-> " << old_time+dt;
      std::cout << ", CFL number = " << cfl << std::endl;

      update_viscosities(dt);

      // update old_solution to current_solution for next time step;
      // this is not done at the end of the previous time step because
      // of the time derivative term in update_viscosities() above
      old_solution = current_solution;

      // solve here
      /** First, compute each \f$\mathbf{f}_i\f$: */
      for (int i = 0; i < rk.s; ++i)
      {
         // compute stage time
         current_time = old_time + rk.c[i]*dt;

         /** compute intermediate solution \f$\mathbf{Y}_i\f$: */
         system_rhs = 0.0;
         mass_matrix.vmult(system_rhs, old_solution);
         for (int j = 0; j < i; ++j)
            system_rhs.add(dt * rk.a[i][j] , rk.f[j]);
         mass_matrix_solve(current_solution);
         // ordinarily, Dirichlet BC need not be reapplied, but in this case, the Dirichlet BC can be time-dependent
         apply_Dirichlet_BC(current_time);

         // residual from solution of previous step is reused in first stage (unless this is the first time step)
         if ((n == 1)||(i != 0))
            compute_ss_residual(rk.f[i]);
      }
      /** Now we compute the solution using the computed \f$\mathbf{f}_i\f$:
       *  \f[
       *    \mathbf{M}\mathbf{y}_{n+1}=\mathbf{M}\mathbf{y}_n + h\sum\limits^s_{i=1}b_i \mathbf{f}_i
       *  \f]
       */    
      if (rk.solution_computed_in_last_stage)
      {
         // end of step residual can be used as beginning of step residual for next step
         rk.f[0] = rk.f[rk.s-1];
      }
      else
      {
         system_rhs = 0.0;
         mass_matrix.vmult(system_rhs, old_solution);
         for (int i = 0; i < rk.s; ++i)
            system_rhs.add(dt * rk.b[i], rk.f[i]);
         mass_matrix_solve(current_solution);
         apply_Dirichlet_BC(current_time);

         // end of step residual can be used as beginning of step residual for next step
         compute_ss_residual(rk.f[0]);
      }
   
      // increment time
      old_time += dt;
      n++;
   
      // output solution of this time step if user has specified;
      //  non-positive numbers for the output_period parameter specify
      //  that solution is not to be output; only the final solution
      //  will be output
      if (conservation_law_parameters.output_period > 0)
      {
         if (n >= next_time_step_output)
         {
            output_solution ();
            next_time_step_output += conservation_law_parameters.output_period;
         }
      }
      else
         if (!(final_time_not_reached))
            output_solution();

      check_nan();

      // compute error for adaptive mesh refinement
      compute_error_for_refinement();

   }// end of time loop
}

/** \fn double ConservationLaw<dim>::compute_dt_from_cfl_condition()
 *  \brief Computes time step size using the CFL condition
 *
 *  The CFL condition for stability is the following:
 *  \f[
 *    \nu = \left|\frac{\lambda_{max}\Delta t}{\Delta x_{min}}\right|\le 1,
 *  \f]
 *  where \f$\lambda_{max}\f$ is the maximum speed in the domain,
 *  \f$\Delta t\f$ is the time step size, and \f$\Delta x_{min}\f$
 *  is the minimum mesh size in the domain. The user supplies the
 *  CFL number (which must be less than 1), and the time step
 *  size is calculated from the CFL definition above.
 */
template <int dim>
double ConservationLaw<dim>::compute_dt_from_cfl_condition()
{
   return conservation_law_parameters.cfl * minimum_cell_diameter / max_flux_speed;
}

/** \fn double ConservationLaw<dim>::compute_cfl_number(const double &dt) const
 *  \brief Computes the CFL number.
 *
 *  The CFL number is the following:
 *  \f[
 *    \nu = \left|\frac{\lambda_{max}\Delta t}{\Delta x_{min}}\right|,
 *  \f]
 *  where \f$\lambda_{max}\f$ is the maximum speed in the domain,
 *  \f$\Delta t\f$ is the time step size, and \f$\Delta x_{min}\f$
 *  is the minimum mesh size in the domain.
 *  \param dt time step size
 *  \return CFL number
 */
template <int dim>
double ConservationLaw<dim>::compute_cfl_number(const double &dt) const
{
   return dt * max_flux_speed / minimum_cell_diameter;
}

/** \fn ConservationLaw<dim>::compute_ss_residual(Vector<double> &f)
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
 *  \param f steady-state residual
 */
template <int dim>
void ConservationLaw<dim>::compute_ss_residual(Vector<double> &f)
{
   // reset vector
   f = 0.0;

   FEValues<dim>     fe_values      (fe, cell_quadrature,
                            update_values | update_gradients | update_JxW_values);
   FEFaceValues<dim> fe_face_values (fe, face_quadrature,
                            update_values | update_gradients | update_JxW_values | update_normal_vectors);

   // allocate memory needed for cell residual and aggregation into global residual
   Vector<double> cell_residual(dofs_per_cell);
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      // reset cell residual
      cell_residual = 0;

      // compute cell contribution to cell residual
      compute_cell_ss_residual(fe_values,
                               cell,
                               cell_residual);

      // compute face contribution to face residual
      /*
      if (need_to_compute_face_residual)
         compute_face_ss_residual(fe_face_values,
                                  cell,
                                  cell_residual);
*/

      // aggregate local residual into global residual
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_residual, local_dof_indices, f);
   } // end cell loop
}

/** \fn void ConservationLaw<dim>::mass_matrix_solve(Vector<double> &x)
 *  \brief Inverts the mass matrix implicitly.
 *
 *  This function computes the product \f$M^{-1}b\f$ of the inverse of the
 *  mass matrix and a vector by solving the linear system \f$M x = b\f$.
 *  The method of inverting the mass matrix is determined by user input.
 *  \param b vector to which the inverse mass matrix is applied
 *  \param x the product \f$M^{-1}b\f$
 */
template <int dim>
void ConservationLaw<dim>::mass_matrix_solve(Vector<double> &x)
{
   linear_solve(conservation_law_parameters.mass_matrix_linear_solver,
                mass_matrix,
                system_rhs,
                x);
}

/** \fn void ConservationLaw<dim>::linear_solve(const typename ConservationLawParameters<dim>::LinearSolverType &linear_solver,
 *                                              const SparseMatrix<double> &A,
 *                                              const Vector<double>       &b,
 *                                                    Vector<double>       &x)
 *  \brief Solves the linear system \f$A x = b\f$.
 *  \param linear_solver linear solution technique to be used to solve the system
 *  \param A the system matrix
 *  \param b the right-hand side vector
 *  \param x the solution vector
 */
template <int dim>
void ConservationLaw<dim>::linear_solve (const typename ConservationLawParameters<dim>::LinearSolverType &linear_solver,
                                         const SparseMatrix<double> &A,
                                         const Vector<double>       &b,
                                               Vector<double>       &x)
{
   switch (linear_solver)
   {
      case ConservationLawParameters<dim>::direct:
      {
         SparseDirectUMFPACK A_umfpack;
         A_umfpack.initialize(A);
         A_umfpack.vmult(x,b);
         break;
      }
      case ConservationLawParameters<dim>::gmres:
      {
         Assert(false,ExcNotImplemented());
         break;
      }
      case ConservationLawParameters<dim>::cg:
      {
         SolverControl solver_control (conservation_law_parameters.max_linear_iterations,
                                       conservation_law_parameters.linear_atol);
         SolverCG<> solver(solver_control);
         
         PreconditionSSOR<> preconditioner;
         preconditioner.initialize(A,1.2); // the second parameter is a relaxation parameter between 1 and 2

         solver.solve(A,x,b,preconditioner);
         break;
      }
      default:
      {
         // throw exception if case was not found
         Assert (false, ExcNotImplemented());
         break;
      }
   }

   // distribute constraints
   constraints.distribute(x);
}

/** \fn void ConservationLaw<dim>::update_viscosities(const double &dt)
 *  \brief Updates viscosity at each quadrature point in each cell.
 */
template <int dim>
void ConservationLaw<dim>::update_viscosities(const double &dt)
{
   switch (conservation_law_parameters.viscosity_type)
   {
      // no viscosity
      case ConservationLawParameters<dim>::none:
      {
         Vector<double> const_visc(n_q_points_cell);
         typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                        endc = dof_handler.end();
         for (; cell != endc; ++cell)
            viscosity_cell_q[cell] = const_visc;
         break;
      }
      // constant viscosity
      case ConservationLawParameters<dim>::constant:
      {
         Vector<double> const_visc(n_q_points_cell);
         const_visc.add(conservation_law_parameters.constant_viscosity_value);
         typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                        endc = dof_handler.end();
         for (; cell != endc; ++cell)
            viscosity_cell_q[cell] = const_visc;
         break;
      }
      // first order viscosity
      case ConservationLawParameters<dim>::first_order:
      {
         update_first_order_viscosities();
         viscosity_cell_q = first_order_viscosity_cell_q;
         break;
      }
      // entropy viscosity
      case ConservationLawParameters<dim>::entropy:
      {
         update_first_order_viscosities();
         update_entropy_viscosities(dt);

         typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                        endc = dof_handler.end();
         if (conservation_law_parameters.add_jumps)
            for (; cell != endc; ++cell)
               for (unsigned int q = 0; q < n_q_points_cell; ++q)
                  viscosity_cell_q[cell](q) = std::min(first_order_viscosity_cell_q[cell](q),
                                                       entropy_viscosity_with_jumps_cell_q[cell](q));
         else
            for (; cell != endc; ++cell)
               for (unsigned int q = 0; q < n_q_points_cell; ++q)
                  viscosity_cell_q[cell](q) = std::min(first_order_viscosity_cell_q[cell](q),
                                                       entropy_viscosity_cell_q[cell](q));
         break;
      }
      default:
      {
         Assert(false,ExcNotImplemented());
         break;
      }
   }
}

/** \fn void ConservationLaw<dim>::update_first_order_viscosities()
 *  \brief Computes first order viscosity at each quadrature point in each cell.
 */
template <int dim>
void ConservationLaw<dim>::update_first_order_viscosities()
{
   double c_max = conservation_law_parameters.first_order_viscosity_coef;

   // loop over cells to compute first order viscosity at each quadrature point
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      double aux = std::abs(c_max * cell_diameter[cell] * max_flux_speed_cell[cell]);
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         first_order_viscosity_cell_q[cell](q) = aux;
   }
}

/** \fn void ConservationLaw<dim>::update_entropy_viscosities(const double &dt)
 *  \brief Computes entropy viscosity at each quadrature point in each cell.
 */
template <int dim>
void ConservationLaw<dim>::update_entropy_viscosities(const double &dt)
{
   // update entropy residuals and max entropy deviation
   update_entropy_residuals(dt);

   // compute entropy viscosity
   double c_s = conservation_law_parameters.entropy_viscosity_coef;
   double c_j = conservation_law_parameters.jump_coef;

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   if (conservation_law_parameters.add_jumps)
   {
      // update jumps
      update_jumps();

      for (; cell != endc; ++cell)
         for (unsigned int q = 0; q < n_q_points_cell; ++q)
         {
            entropy_viscosity_with_jumps_cell_q[cell](q) = std::pow(cell_diameter[cell],2)
               * (c_s * max_entropy_residual_cell[cell] + c_j * max_jumps_cell[cell])
               / max_entropy_deviation;
            // compute entropy viscosity without jumps as well for plotting
            entropy_viscosity_cell_q[cell](q) = std::pow(cell_diameter[cell],2)
               * c_s * max_entropy_residual_cell[cell]
               / max_entropy_deviation;
         }
   } else {
      for (; cell != endc; ++cell)
         for (unsigned int q = 0; q < n_q_points_cell; ++q)
            entropy_viscosity_cell_q[cell](q) = std::pow(cell_diameter[cell],2)
               * c_s * max_entropy_residual_cell[cell]
               / max_entropy_deviation;
   }
}

/** \fn void ConservationLaw<dim>::update_entropy_residuals(const double &dt)
 *  \brief Updates the entropy residuals at each quadrature point in each cell.
 */
template <int dim>
void ConservationLaw<dim>::update_entropy_residuals(const double &dt)
{
   FEValues<dim> fe_values (fe, cell_quadrature, update_values | update_gradients | update_JxW_values);

   std::vector<double>         current_solution_local  (n_q_points_cell);
   std::vector<Tensor<1,dim> > current_gradient_local  (n_q_points_cell);
   Vector<double>              divergence_entropy_flux (n_q_points_cell);

   std::vector<double>         old_solution_local     (n_q_points_cell);
   std::vector<Tensor<1,dim> > old_gradient_local     (n_q_points_cell);
   Vector<double>              old_entropy            (n_q_points_cell);

   // domain-averaged entropy
   double entropy_average = 0.0;

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      fe_values.reinit(cell);
      fe_values.get_function_values(    current_solution, current_solution_local);
      fe_values.get_function_values(    old_solution,     old_solution_local);
      fe_values.get_function_gradients( current_solution, current_gradient_local);
      fe_values.get_function_gradients( old_solution,     old_gradient_local);
      // compute entropy of current and old solutions
      compute_entropy                 (current_solution, fe_values, entropy_cell_q[cell]);
      compute_divergence_entropy_flux (current_solution, fe_values, divergence_entropy_flux);

      compute_entropy            (old_solution,     fe_values, old_entropy);

      for (unsigned int q = 0; q < n_q_points_cell; ++q)
      {
         // compute entropy residual
         double dsdt = (entropy_cell_q[cell](q) - old_entropy(q)) / dt;
         entropy_residual_cell_q[cell](q) = std::abs( dsdt + divergence_entropy_flux(q) );

         // add entropy to volume-weighted sum for use in computation of entropy average
         entropy_average += entropy_cell_q[cell](q) * fe_values.JxW(q);
      }

      // get the max entropy residual on each cell
      max_entropy_residual_cell[cell] = *std::max_element(entropy_residual_cell_q[cell].begin(),entropy_residual_cell_q[cell].end());
   }
   // finish computing entropy average by dividing by the domain volume
   entropy_average /= domain_volume;

   // compute max entropy deviation to be used as entropy normalization term
   max_entropy_deviation = 0.0;
   for (cell = dof_handler.begin_active(); cell != endc; ++cell)
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         max_entropy_deviation = std::max(max_entropy_deviation,
                                          std::abs(entropy_cell_q[cell](q) - entropy_average));
}

/** \fn void ConservationLaw<dim>::update_jumps()
 *  \brief Update the jumps.
 */
template <int dim>
void ConservationLaw<dim>::update_jumps()
{
   FEFaceValues<dim> fe_values_face         (fe, face_quadrature, update_values | update_gradients | update_JxW_values | update_normal_vectors);
   FEFaceValues<dim> fe_values_face_neighbor(fe, face_quadrature, update_values | update_gradients | update_JxW_values | update_normal_vectors);

   std::vector<Tensor<1,dim> > gradients_face         (n_q_points_face);
   std::vector<Tensor<1,dim> > gradients_face_neighbor(n_q_points_face);
   std::vector<Point<dim> >    normal_vectors         (n_q_points_face);
   Vector<double> entropy(n_q_points_face);

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      double max_jump_in_cell = 0.0;
      double max_jump_on_face = 0.0;

      for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
      {
         typename DoFHandler<dim>::face_iterator face = cell->face(iface);
         if (face->at_boundary() == false)
         {
            Assert(cell->neighbor(iface).state() == IteratorState::valid, ExcInternalError());
            typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(iface);
            const unsigned int ineighbor = cell->neighbor_of_neighbor(iface);
            Assert(ineighbor < faces_per_cell, ExcInternalError());

            fe_values_face.reinit(         cell,    iface);
            fe_values_face_neighbor.reinit(neighbor,ineighbor);

            // compute entropy at each quadrature point on face
            compute_entropy_face(current_solution,fe_values_face,entropy);

            // get gradients on adjacent faces of current cell and neighboring cell
            fe_values_face.get_function_gradients(         current_solution, gradients_face);
            fe_values_face_neighbor.get_function_gradients(current_solution, gradients_face_neighbor);

            // get normal vectors
            normal_vectors = fe_values_face.get_normal_vectors();

            max_jump_on_face = 0.0;
            for (unsigned int q = 0; q < n_q_points_face; ++q)
            {
               // compute difference in gradients across face
               gradients_face[q] -= gradients_face_neighbor[q];
               double jump_on_face = std::abs((gradients_face[q] * normal_vectors[q])*entropy[q]);
               max_jump_on_face = std::max(max_jump_on_face, jump_on_face);
            }
         } // end if (at_boundary())
         max_jump_in_cell = std::max(max_jump_in_cell, max_jump_on_face);
      } // end face loop

      max_jumps_cell[cell] = max_jump_in_cell;

   } // end cell loop
}

/** \fn void ConservationLaw<dim>::check_nan()
 *  \brief Checks that there are no NaNs in the solution vector
 *
 *  The NaN check is performed by comparing a value to itself
 *  since an equality comparison with NaN always returns false.
 */
template <int dim>
void ConservationLaw<dim>::check_nan()
{
   unsigned int n = dof_handler.n_dofs();
   for (unsigned int i = 0; i < n; ++i)
      Assert(current_solution(i) == current_solution(i), ExcNumberNotFinite());
}
