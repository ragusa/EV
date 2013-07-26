/** \fn ConservationLaw
 *  \brief Constructor for ConservationLaw class.
 * 
 *  In addition to initializing some member variables,
 *  this funciton gets parameters from the parameter
 *  handler.
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

/** \fn run
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

/** \fn invert_mass_matrix
 *  \brief Inverts the mass matrix implicitly.
 *
 *  This function implicitly inverts the mass matrix by solving
 *  the linear system M*x=b
 */
template <int dim>
Vector<double> ConservationLaw<dim>::invert_mass_matrix(Vector<double> b)
{
   // perform direct solve using UMFPACK
   Vector<double> result;
   return result;
}

/** \fn solve_erk
 *  \brief Solves the transient using explicit Runge-Kutta.
 *
 *  This function contains the transient loop and solves the
 *  transient using explicit Runge-Kutta
 */
template <int dim>
void ConservationLaw<dim>::solve_erk()
{
   // get ERK parameters
   int Ns = conservation_law_parameters.erk_nstages;
   std::vector<Vector<double> > erk_a(Ns);
   for (int i = 0; i < Ns; ++i)
      erk_a[i].reinit(Ns);
   std::vector<double> erk_b(Ns);
   std::vector<double> erk_c(Ns);
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
      erk_k[0] = compute_ss_residual(time,old_solution);
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

         erk_k[i] = compute_ss_residual(time + erk_c[i]*dt, y_tmp);
      }
      // compute new solution
      current_solution = old_solution;
      for (int i = 0; i < Ns; ++i)
      {
         x_tmp = erk_k[i];
         x_tmp *= (dt * erk_b[i]);
         current_solution += x_tmp;
      }
   
      // increment time
      time += dt;
   
      // output solution of this time step if user has specified;
      //  negative numbers for the output_period parameter specify
      //  that solution is to be output every time step
      if (conservation_law_parameters.output_period < 0)
         output_results ();
      else if (time >= next_time_step_output)
      {
         output_results ();
         next_time_step_output += conservation_law_parameters.output_period;
      }

      // update old_solution to current_solution for next time step
      old_solution = current_solution;
   }
}

/** \fn setup_system
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
   right_hand_side.reinit(dof_handler.n_dofs());

}

/** \fn compute_ss_residual
 *  \brief Computes the steady state residual.
 *
 *  This function is to be defined in the derived physics class.
 */
template <int dim>
Vector<double> ConservationLaw<dim>::compute_ss_residual (double t, Vector<double> solution)
{
   Vector<double> some_vector(dof_handler.n_dofs());
   return some_vector;
}

/*
template <int dim>
void ConservationLaw<dim>::assemble_cell_term (const FEValues<dim>             &fe_v,
                                               const std::vector<unsigned int> &dofs)
{}

template <int dim>
void ConservationLaw<dim>::assemble_face_term (const unsigned int               face_no,
                                               const FEFaceValuesBase<dim>     &fe_v,
                                               const FEFaceValuesBase<dim>     &fe_v_neighbor,
                                               const std::vector<unsigned int> &dofs,
                                               const std::vector<unsigned int> &dofs_neighbor,
                                               const bool                       external_face,
                                               const unsigned int               boundary_id,
                                               const double                     face_diameter)
{}
*/

/** \fn linear_solve
 *  \brief Performs a linear solve.
 *
 *  This is only called if the time integrator is implicit. The
 *  function contains a switch for the linear solver type.
 */
/*
template <int dim>
std::pair<unsigned int, double> ConservationLaw<dim>::linear_solve (Vector<double> &newton_update)
{
	switch (conservation_law_parameters.linear_solver)
	{
	case ConservationLawParameters<dim>::direct:
				{
		//SparseDirectUMFPACK A_direct;
		//A_direct.initialize(system_matrix);
		//A_direct.vmult(newton_update, right_hand_side);

		Assert(false,ExcNotImplemented());
		return std::pair<unsigned int, double> (0,0);
				}
	case ConservationLawParameters<dim>::gmres:
	{
            Assert(false,ExcNotImplemented());
            return std::pair<unsigned int, double> (0,0);
	}
	case ConservationLawParameters<dim>::bicgstab:
	{
        SolverControl solver_control(1000, 1e-6);
        SolverBicgstab<> solver(solver_control);

        solver.solve(system_matrix, newton_update, right_hand_side,
                    PreconditionIdentity());
        return std::pair<unsigned int, double> (solver_control.last_step(),
        		solver_control.last_value());
            //Assert(false,ExcNotImplemented());
            //return std::pair<unsigned int, double> (0,0);
	}
	}

	// throw exception if case was not found
    Assert (false, ExcNotImplemented());
    return std::pair<unsigned int, double> (0,0);
}
*/

/*
template <int dim>
void ConservationLaw<dim>::compute_refinement_indicators (Vector<double> &indicator) const
{}

template <int dim>
void ConservationLaw<dim>::refine_grid (const Vector<double> &indicator)
{}
*/

/** \fn output_results
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
	std::string filename = "solution-" +
			Utilities::int_to_string (output_file_number, 3) +
			".vtk";
	std::ofstream output (filename.c_str());
	data_out.write_vtk (output);

	++output_file_number;
}
