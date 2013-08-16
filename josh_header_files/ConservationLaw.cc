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
   quadrature(2*params.degree + 1),
   n_q_points_cell(quadrature.size()),
   face_quadrature(2*params.degree),
   n_q_points_face(face_quadrature.size()),
   verbose_cout(std::cout, false),
   initial_conditions(params.n_components)
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
   // setup system
   setup_system();

   // interpolate the initial conditions to the grid
   VectorTools::interpolate(dof_handler,initial_conditions,current_solution);

   // apply Dirichlet BC to initial solution or guess
   apply_Dirichlet_BC();
   
   old_solution = current_solution;

   // output initial solution
   output_solution();

   // begin time stepping loop
   switch (conservation_law_parameters.temporal_integrator)
   {
      case ConservationLawParameters<dim>::erk: // explicit Runge-Kutta
          solve_erk();
          break;
      default:
          Assert(false,ExcNotImplemented());
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

/** \fn void ConservationLaw<dim>::setup_system()
 *  \brief Sets up the system before solving.
 *
 *  This function makes the sparsity pattern and reinitializes
 *  the system matrix with the sparsity pattern.
 */
template <int dim>
void ConservationLaw<dim>::setup_system ()
{
   // clear maps
   dx.clear();
   flux_speed_cell_q.clear();
   max_flux_speed_cell.clear();
   viscosity_cell_q.clear();
   first_order_viscosity_cell_q.clear();
   entropy_viscosity_cell_q.clear();
   entropy_viscosity_with_jumps_cell_q.clear();
   entropy_residual_cell_q.clear();
   max_jumps_cell.clear();

   // make grid and refine
   double domain_start = 0;
   double domain_width = conservation_law_parameters.domain_width;
   domain_volume = std::pow(domain_width,dim);
   GridGenerator::hyper_cube(triangulation, domain_start, domain_start + domain_width);
   triangulation.refine_global(conservation_law_parameters.initial_refinement_level);
   update_cell_sizes();

   // clear and distribute dofs
   dof_handler.clear();
   dof_handler.distribute_dofs(fe);

   // make constraints
   constraints.clear();
   DoFTools::make_hanging_node_constraints (dof_handler, constraints);
   for (int boundary = 0; boundary < conservation_law_parameters.n_boundaries; ++boundary)
   {
      for (int component = 0; component < n_components; ++component)
      {
         bool is_Dirichlet = true;
         if (is_Dirichlet)
         {
            std::vector<bool> component_mask(n_components, false);
            component_mask[component] = true;
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      boundary,
                                                      ZeroFunction<dim>(n_components),
                                                      constraints,
                                                      component_mask);
         }
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

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      flux_speed_cell_q[cell]            = Vector<double>(n_q_points_cell);
      viscosity_cell_q[cell]             = Vector<double>(n_q_points_cell);
      first_order_viscosity_cell_q[cell] = Vector<double>(n_q_points_cell);
      entropy_viscosity_cell_q[cell]     = Vector<double>(n_q_points_cell);
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
   // reset minimum cell size
   dx_min = conservation_law_parameters.domain_width;

   // fill cell size map and find minimum cell size
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      dx[cell] = cell->diameter();
      dx_min = std::min( dx_min, dx[cell] );
   }
}

/** \fn void ConservationLaw<dim>::assemble_mass_matrix()
 *  \brief Assembles the mass matrix and applies constraints.
 */
template <int dim>
void ConservationLaw<dim>::assemble_mass_matrix ()
{
   mass_matrix = 0.0;

   FEValues<dim> fe_values (fe, quadrature, update_values | update_gradients | update_JxW_values);

   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      FullMatrix<double> local_mass (dofs_per_cell, dofs_per_cell);

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

   // output mass matrix
   if (conservation_law_parameters.output_mass_matrix)
   {
      std::ofstream mass_matrix_out ("mass_matrix.txt");
      mass_matrix.print_formatted(mass_matrix_out, 10, true, 0, "0", 1);
      mass_matrix_out.close();
   }
}

/** \fn void ConservationLaw<dim>apply_Dirichlet_BC()
 *  \brief Applies Dirichlet boundary conditions
 *
 *  This function applies Dirichlet boundary conditions
 *  using the interpolate_boundary_values() tool.
 */
template <int dim>
void ConservationLaw<dim>::apply_Dirichlet_BC()
{
   std::map<unsigned int, double> boundary_values;
   for (int boundary = 0; boundary < conservation_law_parameters.n_boundaries; ++boundary)
      for (int component = 0; component < n_components; ++component)
      {
         bool is_dirichlet = true;
         if (is_dirichlet)
         {
            std::vector<bool> component_mask(n_components, false);
            component_mask[component] = true;
            VectorTools::interpolate_boundary_values (dof_handler,
                                                      boundary,
                                                      ZeroFunction<dim>(),
                                                      boundary_values,
                                                      component_mask);
         }
      }
   for (std::map<unsigned int, double>::const_iterator it = boundary_values.begin(); it != boundary_values.end(); ++it)
      current_solution(it->first) = (it->second);
}

/** \fn void ConservationLaw<dim>::output_solution() const
 *  \brief Outputs the solution to .vtk.
 *
 *  The user supplies an input parameter that determines
 *  how often this function is called in the transient.
 *  A static variable for the output file number is incremented
 *  the function is called.
 */
template <int dim>
void ConservationLaw<dim>::output_solution () const
{
   DataOut<dim> data_out;
   data_out.attach_dof_handler (dof_handler);

   data_out.add_data_vector (current_solution,
                             component_names,
                             DataOut<dim>::type_dof_data,
                             component_interpretations);

   data_out.add_data_vector (current_solution, "current_solution");

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

/** \fn void ConservationLaw<dim>::output_map(const std::map<DoFHandler<dim>::active_cell_iterator &map,
                                              const std::string &output_filename_base)
 *  \brief Outputs a mapped quantity at all quadrature points
 */
template <int dim>
void ConservationLaw<dim>::output_map(std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > &map,
                                      const std::string &output_filename_base)
{
   if (dim == 1)
   {
      // get date from maps into vector of point-data pairs
      unsigned int n_cells = triangulation.n_active_cells();
      unsigned int total_n_q_points = n_cells * n_q_points_cell;
      std::vector<std::pair<double,double> > profile(total_n_q_points);

      FEValues<dim> fe_values (fe, quadrature, update_quadrature_points);
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
      Assert(false, ExcNotImplemented());
   }
}

/** \fn ConservationLaw<dim>::solve_erk()
 *  \brief Solves the transient using explicit Runge-Kutta.
 *
 *  This function contains the transient loop and solves the
 *  transient using explicit Runge-Kutta:
 *  \f[
 *    y_{n+1}=y_n + \sum\limits^s_{i=1}b_i M^{-1} k_i
 *  \f]
 *  where \f$k_i\f$ is
 *  \f[
 *    k_i=h f(t_n + c_i h, y_n + \sum\limits^{i-1}_{j=1}a_{i,j} M^{-1} k_j)
 *  \f]
 */
template <int dim>
void ConservationLaw<dim>::solve_erk()
{
   // get ERK parameters a, b, and c (Butcher tableau)
   int Ns = conservation_law_parameters.erk_nstages;
   std::vector<Vector<double> > erk_a(Ns);
   for (int i = 0; i < Ns; ++i)
      erk_a[i].reinit(Ns);
   std::vector<double> erk_b(Ns);
   std::vector<double> erk_c(Ns);
   // for now, assume there is only one ERK method for each number of stages
   switch (Ns)
   {
      case 1:
         erk_b[0] = 1;
         erk_c[0] = 0;
         break;
      case 2:
         erk_a[1][0] = 0.5;
         erk_b[0] = 0;
         erk_b[1] = 1;
         erk_c[0] = 0;
         erk_c[1] = 0.5;
         break;
      case 3:
         erk_a[1][0] = 1.0;
         erk_a[2][0] = 0.25;
         erk_a[2][1] = 0.25;
         erk_b[0] = 1./6;
         erk_b[1] = 1./6;
         erk_b[2] = 4./6;
         erk_c[0] = 0;
         erk_c[1] = 1.0;
         erk_c[2] = 0.5;
         break;
      case 4:
         erk_a[1][0] = 0.5;
         erk_a[2][0] = 0;
         erk_a[2][1] = 0.5;
         erk_a[3][0] = 0;
         erk_a[3][1] = 0;
         erk_a[3][2] = 1;
         erk_b[0] = 1./6;
         erk_b[1] = 1./3;
         erk_b[2] = 1./3;
         erk_b[3] = 1./6;
         erk_c[0] = 0;
         erk_c[1] = 0.5;
         erk_c[2] = 0.5;
         erk_c[3] = 1;
         break;
      default:
         Assert(false,ExcNotImplemented());
         break;
   }

   // allocate memory for each stage
   std::vector<Vector<double> > erk_f(Ns);
   for (int i = 0; i < Ns; ++i)
      erk_f[i].reinit(dof_handler.n_dofs());

   old_time = 0.0;
   unsigned int n = 1; // time step index
   unsigned int next_time_step_output = conservation_law_parameters.output_period;
   double t_end = conservation_law_parameters.final_time;
   bool final_time_not_reached = true;
   while (final_time_not_reached)
   {
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

      std::cout << "time step " << n << ": t = " << old_time << "-> " << old_time+dt << std::endl;
      std::cout << "   Number of active cells: ";
      std::cout << triangulation.n_active_cells();
      std::cout << std::endl;
      std::cout << "   Number of degrees of freedom: ";
      std::cout << dof_handler.n_dofs();
      std::cout << std::endl;
      std::cout << std::endl;

      // solve here
      /** First, compute each \f$\mathbf{M}^{-1}\mathbf{k}_i\f$: */
      for (int i = 0; i < Ns; ++i)
      {
         // compute stage time
         current_time = old_time + erk_c[i]*dt;

         /** compute intermediate solution: \f$\mathbf{y}_n + \sum\limits^{i-1}_{j=1}a_{i,j} \mathbf{M}^{-1} \mathbf{k}_j\f$ */
         system_rhs = 0.0;
         mass_matrix.vmult(system_rhs, old_solution);
         for (int j = 0; j < i-1; ++j)
            system_rhs.add(dt * erk_a[i][j] , erk_f[j]);
         mass_matrix_solve(current_solution);
         apply_Dirichlet_BC();

         compute_ss_residual(erk_f[i]);
      }
      /** Now we compute the solution using the computed \f$\mathbf{k}_i\f$:
       *  \f[
       *    \mathbf{y}_{n+1}=\mathbf{y}_n + \sum\limits^s_{i=1}b_i \mathbf{M}^{-1} \mathbf{k}_i
       *  \f]
       */    
      system_rhs = 0.0;
      mass_matrix.vmult(system_rhs, old_solution);
      for (int i = 0; i < Ns; ++i)
         system_rhs.add(dt * erk_b[i], erk_f[i]);
      mass_matrix_solve(current_solution);
      apply_Dirichlet_BC();
   
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

      // update old_solution to current_solution for next time step
      old_solution = current_solution;
      check_nan();
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
   // update max speed
   update_flux_speeds();

   // compute time step size from CFL condition
   double dt = conservation_law_parameters.cfl * dx_min / max_flux_speed;

   return dt;
}

/** \fn ConservationLaw<dim>::update_flux_speeds()
 *  \brief Updates the map of flux speeds and recomputes max flux speed.
 */
template <int dim>
void ConservationLaw<dim>::update_flux_speeds()
{
   FEValues<dim> fe_values(fe, quadrature, update_values);
   Tensor<1,dim> dfdu;
   std::vector<double> local_solution(n_q_points_cell);

   // reset max flux speed
   max_flux_speed = 0.0;

   // loop over cells to compute first order viscosity at each quadrature point
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      fe_values.reinit(cell);
      fe_values.get_function_values(current_solution, local_solution);

      for (unsigned int q = 0; q < n_q_points_cell; ++q)
      {
         dfdu = flux_derivative(local_solution[q]);
         flux_speed_cell_q[cell](q) = dfdu.norm();
      }

      // get max flux speed
      max_flux_speed_cell[cell] = *std::max_element(flux_speed_cell_q[cell].begin(), flux_speed_cell_q[cell].end());
      max_flux_speed = std::max(max_flux_speed, max_flux_speed_cell[cell]);
   }
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

   // update viscosities
   update_viscosities();

   FEValues<dim>     fe_values      (fe, quadrature,
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
      compute_face_ss_residual(fe_face_values,
                               cell,
                               cell_residual);

      // aggregate local residual into global residual
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_residual, local_dof_indices, f);
   } // end cell loop
}

/** \fn void ConservationLaw<dim>::mass_matrix_solve(const Vector<double> &b,
 *                                                         Vector<double> &x)
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

/** \fn void ConservationLaw<dim>::update_viscosities(const double dt)
 *  \brief Updates viscosity at each quadrature point in each cell.
 */
template <int dim>
void ConservationLaw<dim>::update_viscosities()
{
   // update flux speeds
   update_flux_speeds();

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
         update_entropy_viscosities();

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
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         first_order_viscosity_cell_q[cell](q) = std::abs(c_max * dx[cell] * flux_speed_cell_q[cell](q));
}

/** \fn void ConservationLaw<dim>::update_entropy_viscosities()
 *  \brief Computes entropy viscosity at each quadrature point in each cell.
 */
template <int dim>
void ConservationLaw<dim>::update_entropy_viscosities()
{
   // update entropy residuals and max entropy deviation
   update_entropy_residuals();

   // compute entropy viscosity
   double c_s = conservation_law_parameters.entropy_viscosity_coef;

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   if (conservation_law_parameters.add_jumps)
   {
      // update jumps
      update_jumps();

      for (; cell != endc; ++cell)
         for (unsigned int q = 0; q < n_q_points_cell; ++q)
         {
            entropy_viscosity_with_jumps_cell_q[cell](q) = c_s * std::pow(dx[cell],2)
               * (entropy_residual_cell_q[cell](q) + max_jumps_cell[cell])
               / max_entropy_deviation;
            // compute entropy viscosity without jumps as well for plotting
            entropy_viscosity_cell_q[cell](q) = c_s * std::pow(dx[cell],2)
               * entropy_residual_cell_q[cell](q)
               / max_entropy_deviation;
         }
   } else {
      for (; cell != endc; ++cell)
         for (unsigned int q = 0; q < n_q_points_cell; ++q)
            entropy_viscosity_cell_q[cell](q) = c_s * std::pow(dx[cell],2)
               * entropy_residual_cell_q[cell](q)
               / max_entropy_deviation;
   }
}

/** \fn void ConservationLaw<dim>::update_entropy_residuals(const double dt)
 *  \brief Updates the entropy residuals at each quadrature point in each cell.
 */
template <int dim>
void ConservationLaw<dim>::update_entropy_residuals()
{
   FEValues<dim> fe_values (fe, quadrature, update_values | update_gradients | update_JxW_values);

   std::vector<double> current_solution_local(n_q_points_cell);
   std::vector<double> old_solution_local    (n_q_points_cell);
   std::vector<Tensor<1,dim> > current_gradient_local(n_q_points_cell);
   std::vector<Tensor<1,dim> > old_gradient_local    (n_q_points_cell);
   std::vector<Tensor<1,dim> > current_flux_derivative(n_q_points_cell);
   std::vector<Tensor<1,dim> > old_flux_derivative    (n_q_points_cell);

   std::vector<double> current_entropy(n_q_points_cell);
   std::vector<double> old_entropy    (n_q_points_cell);
   std::vector<double> current_entropy_derivative(n_q_points_cell);
   std::vector<double> old_entropy_derivative    (n_q_points_cell);

   double dt = current_time - old_time;

   // domain-averaged entropy
   double entropy_average = 0.0;

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      fe_values.reinit(cell);
      fe_values.get_function_values(current_solution, current_solution_local);
      fe_values.get_function_values(old_solution,     old_solution_local);
      fe_values.get_function_gradients(current_solution, current_gradient_local);
      fe_values.get_function_gradients(old_solution,     old_gradient_local);

      for (unsigned int q = 0; q < n_q_points_cell; ++q)
      {
         current_entropy[q] = entropy(current_solution_local[q]);
         old_entropy[q]     = entropy(old_solution_local[q]);
         current_entropy_derivative[q] = entropy_derivative(current_solution_local[q]);
         old_entropy_derivative[q]     = entropy_derivative(old_solution_local[q]);
         current_flux_derivative[q] = flux_derivative(current_solution_local[q]);
         old_flux_derivative[q]     = flux_derivative(old_solution_local[q]);

         // compute entropy residual
         entropy_residual_cell_q[cell](q) = std::abs( (current_solution_local[q] - old_solution_local[q]) / dt
            + current_entropy_derivative[q] * current_flux_derivative[q] * current_gradient_local[q] );

         // add entropy to volume-weighted sum for use in computation of entropy average
         entropy_average += entropy(current_solution_local[q]) * fe_values.JxW(q);
      }
   }
   // finish computing entropy average by dividing by the domain volume
   entropy_average /= domain_volume;

   // compute max entropy deviation to be used as entropy normalization term
   max_entropy_deviation = 0.0;
   for (cell = dof_handler.begin_active(); cell != endc; ++cell)
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         max_entropy_deviation = std::max(max_entropy_deviation,
                                          std::abs(current_entropy[q] - entropy_average));
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

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      double max_jump_in_cell = 0.0;
      double max_jump_on_face;

      for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
      {
         typename DoFHandler<dim>::face_iterator face = cell->face(iface);
         if (face->at_boundary() == false)
         {
            Assert(cell->neighbor(iface).state() == IteratorState::valid, ExcInternalError());
            typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(iface);
            const unsigned int ineighbor = cell->neighbor_of_neighbor(iface);
            Assert(ineighbor < faces_per_cell, ExcInternalError());

            fe_values_face.reinit(cell,iface);
            fe_values_face_neighbor.reinit(cell,ineighbor);

            // get gradients on adjacent faces of current cell and neighboring cell
            fe_values_face.get_function_gradients(current_solution, gradients_face);
            fe_values_face_neighbor.get_function_gradients(current_solution, gradients_face_neighbor);

            // get normal vectors
            normal_vectors = fe_values_face.get_normal_vectors();

            max_jump_on_face = 0.0;
            for (unsigned int q = 0; q < n_q_points_face; ++q)
            {
               // compute difference in gradients across face
               gradients_face[q] -= gradients_face_neighbor[q];
               double jump_on_face = std::abs(gradients_face[q] * normal_vectors[q]);
               max_jump_on_face = std::max(max_jump_on_face, jump_on_face);
            }
         } // end if (at_boundary())
         max_jump_in_cell = std::max(max_jump_in_cell, max_jump_on_face);
      } // end face loop

      max_jumps_cell[cell] = max_jump_in_cell;

   } // end cell loop
}

/*
template <int dim>
void ConservationLaw<dim>::compute_refinement_indicators (Vector<double> &indicator) const
{}

template <int dim>
void ConservationLaw<dim>::refine_grid (const Vector<double> &indicator)
{}
*/

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
