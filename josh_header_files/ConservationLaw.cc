/** \file ConservationLaw.cc
 *  \brief Provides function definitions for the ConservationLaw class.
 */

/** \fn ConservationLaw<dim>::ConservationLaw(ParameterHandler &prm, const int &n_comp)
 *  \brief Constructor for ConservationLaw class.
 * 
 *  In addition to initializing some member variables,
 *  this funciton gets parameters from the parameter
 *  handler.
 *  \param prm parameter handler containing conservation law parameters
 *  \param n_comp number of components in the system
 */
template <int dim>
ConservationLaw<dim>::ConservationLaw(ParameterHandler &prm,
                                      const int &n_comp):
   n_components(n_comp),
   mapping(),
   fe(FE_Q<dim>(1), n_comp),
   dof_handler(triangulation),
   quadrature(2),
   face_quadrature(2),
   verbose_cout(std::cout, false),
   initial_conditions(n_comp)
{
   // get conservation law parameters
   conservation_law_parameters.get_parameters(prm);
}

/** \fn ConservationLaw<dim>::run()
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
   VectorTools::interpolate(dof_handler,initial_conditions,old_solution);
   current_solution = old_solution;
   
   // output initial solution
   output_results();

   // begin time stepping loop
   switch (conservation_law_parameters.temporal_integrator)
   {
      case ConservationLawParameters<dim>::erk: // explicit Runge-Kutta
          solve_erk();
          break;
      default:
          Assert(false,ExcNotImplemented());
   }
}

/** \fn void ConservationLaw<dim>::invert_mass_matrix(const Vector<double> &b,
 *                                                          Vector<double> &x)
 *  \brief Inverts the mass matrix implicitly.
 *
 *  This function computes the product \f$M^{-1}b\f$ of the inverse of the
 *  mass matrix and a vector by solving the linear system \f$M x = b\f$.
 *  The method of inverting the mass matrix is determined by user input.
 *  \param b vector to which the inverse mass matrix is applied
 *  \param x the product \f$M^{-1}b\f$
 */
template <int dim>
void ConservationLaw<dim>::invert_mass_matrix(const Vector<double> &b,
                                                    Vector<double> &x)
{
   linear_solve(conservation_law_parameters.mass_matrix_linear_solver,
                mass_matrix,
                b,
                x);
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
         erk_b[0] = 0;
         erk_c[0] = 1;
         break;
      case 2:
         erk_a[1][0] = 0.5;
         erk_b[0] = 0;
         erk_b[1] = 0.5;
         erk_c[0] = 0;
         erk_c[1] = 1;
         break;
      case 4:
         erk_a[1][0] = 0.5;
         erk_a[2][0] = 0;
         erk_a[2][1] = 0.5;
         erk_a[3][0] = 0;
         erk_a[3][1] = 0;
         erk_a[3][2] = 1;
         erk_b[0] = 0;
         erk_b[1] = 0.5;
         erk_b[2] = 0.5;
         erk_b[3] = 1;
         erk_c[0] = 1./6;
         erk_c[1] = 1./3;
         erk_c[2] = 1./3;
         erk_c[3] = 1./6;
         break;
      default:
         Assert(false,ExcNotImplemented());
         break;
   }

   // allocate memory for each stage
   std::vector<Vector<double> > erk_k(Ns);
   for (int i = 0; i < Ns; ++i)
      erk_k[i].reinit(dof_handler.n_dofs());

   // allocate memory for intermediate steps
   Vector<double> y_tmp(dof_handler.n_dofs());
   Vector<double> x_tmp(dof_handler.n_dofs());

   double time = 0;
   unsigned int next_time_step_output = conservation_law_parameters.output_period;
   double dt = conservation_law_parameters.time_step_size;
   double t_end = conservation_law_parameters.final_time;
   bool final_time_not_reached_yet = true;
   while (final_time_not_reached_yet)
   {
      std::cout << "time: " << time << std::endl;
      std::cout << "   Number of active cells: ";
      std::cout << triangulation.n_active_cells();
      std::cout << std::endl;
      std::cout << "   Number of degrees of freedom: ";
      std::cout << dof_handler.n_dofs();
      std::cout << std::endl;
      std::cout << std::endl;

      // check end of transient and shorten last time step if necessary
      if ((time+dt) >= t_end)
      {
         dt = t_end - time;
         final_time_not_reached_yet = false;
      }

      // solve here
      compute_ss_residual(time,old_solution);
      invert_mass_matrix(ss_residual,erk_k[0]);
      erk_k[0] *= dt;
      for (int i = 1; i < Ns; ++i)
      {
         // compute intermediate solution
         y_tmp = old_solution;
         for (int j = 0; j < i-1; ++j)
         {
            x_tmp = erk_k[j];
            x_tmp *= erk_a[i][j];
            y_tmp += x_tmp;
         }

         compute_ss_residual(time + erk_c[i]*dt, y_tmp);
         invert_mass_matrix(ss_residual,erk_k[i]);
         erk_k[i] *= dt;
      }
      // compute new solution
      current_solution = old_solution;
      for (int i = 0; i < Ns; ++i)
      {
         x_tmp = erk_k[i];
         x_tmp *= erk_b[i];
         current_solution += x_tmp;
      }
   
      // increment time
      time += dt;
   
      // output solution of this time step if user has specified;
      //  non-positive numbers for the output_period parameter specify
      //  that solution is not to be output; only the final solution
      //  will be output
      if (conservation_law_parameters.output_period > 0)
         if (time >= next_time_step_output)
         {
            output_results ();
            next_time_step_output += conservation_law_parameters.output_period;
         }
      else
         if (!(final_time_not_reached_yet))
            output_results();

      // update old_solution to current_solution for next time step
      old_solution = current_solution;
   }// end of time loop
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
   // make grid and refine
   GridGenerator::hyper_cube(triangulation,-1,1);
   triangulation.refine_global(3);

   // clear and distribute dofs
   dof_handler.clear();
   dof_handler.distribute_dofs(fe);

   // create sparsity pattern and compress
   sparsity_pattern.reinit (dof_handler.n_dofs(),
                            dof_handler.n_dofs(),
                            dof_handler.max_couplings_between_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
   sparsity_pattern.compress();

   // create mass matrix
   mass_matrix.reinit(sparsity_pattern);
   MatrixCreator::create_mass_matrix( dof_handler,
                                      QGauss<dim>(3),
                                      mass_matrix);

   // resize vectors
   old_solution.reinit(dof_handler.n_dofs());
   current_solution.reinit(dof_handler.n_dofs());
   ss_residual.reinit(dof_handler.n_dofs());

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
/*
      case ConservationLawParameters<dim>::bicgstab:
      {
         SolverControl solver_control(1000, 1e-6);
         SolverBicgstab<> solver(solver_control);

         solver.solve(system_matrix, newton_update, right_hand_side,
                      PreconditionIdentity());
         return std::pair<unsigned int, double> (solver_control.last_step(),
                                                 solver_control.last_value());
      }
*/
      default:
      {
         // throw exception if case was not found
         Assert (false, ExcNotImplemented());
         break;
      }
   }
}

/*
template <int dim>
void ConservationLaw<dim>::compute_refinement_indicators (Vector<double> &indicator) const
{}

template <int dim>
void ConservationLaw<dim>::refine_grid (const Vector<double> &indicator)
{}
*/

/** \fn void ConservationLaw<dim>::output_results() const
 *  \brief Outputs the solution to .vtk.
 *
 *  The user supplies an input parameter that determines
 *  how often this function is called in the transient.
 *  A static variable for the output file number is incremented
 *  the function is called.
 */
template <int dim>
void ConservationLaw<dim>::output_results () const
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
   std::string filename = "output/solution-" +
         Utilities::int_to_string (output_file_number, 3) +
         ".vtk";
   std::ofstream output (filename.c_str());
   data_out.write_vtk (output);

   ++output_file_number;
}
