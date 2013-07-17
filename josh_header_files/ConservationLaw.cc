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

template <int dim>
void ConservationLaw<dim>::run()
{
   // make grid and refine
   GridGenerator::hyper_cube(triangulation,-1,1);
   triangulation.refine_global(3);

   // clear and distribute dofs
   dof_handler.clear();
   dof_handler.distribute_dofs(fe);

   // resize vectors
   old_old_solution.reinit(dof_handler.n_dofs());
   old_solution.reinit(dof_handler.n_dofs());
   current_solution.reinit(dof_handler.n_dofs());
   right_hand_side.reinit(dof_handler.n_dofs());
   predictor.reinit(dof_handler.n_dofs());

   // setup system
   setup_system();

   // interpolate the initial conditions to the grid
   VectorTools::interpolate(dof_handler,initial_conditions,old_solution);
   
   current_solution = old_solution;
   predictor = old_solution;

   assemble_system();

   // output initial solution
   output_results();

   // initialize Newton update
   Vector<double> newton_update(dof_handler.n_dofs());

   // begin time stepping loop
   double time = 0;
   unsigned int next_time_step_output = conservation_law_parameters.output_period;
   while (time < conservation_law_parameters.final_time)
   {
	   std::cout << "time: " << time << std::endl;
	   std::cout << "   Number of active cells: ";
	   std::cout << triangulation.n_active_cells();
	   std::cout << std::endl;
	   std::cout << "   Number of degrees of freedom: ";
	   std::cout << dof_handler.n_dofs();
	   std::cout << std::endl;
	   std::cout << std::endl;

	   std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl
                 << "   _____________________________________" << std::endl;

	   // reset nonlinear iteration number to zero
	   int nonlin_iter = 0;

	   // set guess for this time step's solution to be predictor from last time step
	   current_solution = predictor;

	   // initialize convergence flag
	   bool not_converged_yet = true;
	   // begin Newton loop for current time step
	   while (not_converged_yet)
	   {
		   // reset system matrix and rhs vector and then recompute
		   system_matrix = 0;
		   right_hand_side = 0;
		   assemble_system ();

		   // compute L2 norm of right hand side
		   const double nonlinear_residual_norm = right_hand_side.l2_norm();
		   // determine if nonlinear solution has converged
		   if (std::fabs(nonlinear_residual_norm) < conservation_law_parameters.nonlinear_atol)
		   // nonlinear solution has converged
		   {
			   // print nonlinear residual norm and report convergence
			   std::printf("   %-16.3e (converged)\n\n", nonlinear_residual_norm);
			   break;
		   }
		   else
		   // nonlinear solution has not yet converged; take another Newton step
		   {
			   // reset Newton update
               newton_update = 0;

               // solve linear system
			   std::pair<unsigned int, double> linear_solve_info = solve (newton_update);

			   // update solution with Newton update
			   current_solution += newton_update;

			   // print nonlinear residual norm, number of linear iterations, and linear residual norm
			   if (conservation_law_parameters.linear_solver !=
					   ConservationLawParameters<dim>::direct)
				   std::printf("   %-16.3e %04d        %-5.2e\n",
						   nonlinear_residual_norm, linear_solve_info.first,
						   linear_solve_info.second);
			   else
				   std::printf("   %-16.3e N/A        N/A\n",
						   nonlinear_residual_norm);
		   }

		   // increment nonlinear iteration number
		   ++nonlin_iter;

		   // throw error if the maximum number of nonlinear iterations has been exceeded
		   AssertThrow (nonlin_iter <= conservation_law_parameters.max_nonlinear_iterations,
				   ExcMessage ("No convergence in nonlinear solver"));
	   }
	   // increment time
	   time += conservation_law_parameters.time_step_size;

	   /* output solution of this time step if user has specified;
	      negative numbers for the output_period parameter specify
	      that solution is to be output every time step
	   */
	   if (conservation_law_parameters.output_period < 0)
		   output_results ();
	   else if (time >= next_time_step_output)
	   {
		   output_results ();
		   next_time_step_output += conservation_law_parameters.output_period;
	   }

	   // predict solution for next time step with u^(n+1) ~= 2*u^(n) - u^(n-1)
	   predictor = current_solution;
	   predictor.sadd (2.0, -1.0, old_solution);

	   // update old_solution to current_solution for next time step
	   old_solution = current_solution;
   }
}

template <int dim>
void ConservationLaw<dim>::setup_system ()
{
   CompressedSparsityPattern compressed_sparsity_pattern (dof_handler.n_dofs(),
                                                          dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);
   sparsity_pattern.copy_from(compressed_sparsity_pattern);
 
   system_matrix.reinit (sparsity_pattern);
}

template <int dim>
void ConservationLaw<dim>::assemble_system ()
{
   std::cout << "Assemble system" << std::endl;
}

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

template <int dim>
std::pair<unsigned int, double> ConservationLaw<dim>::solve (Vector<double> &newton_update)
{
	switch (conservation_law_parameters.linear_solver)
	{
	case ConservationLawParameters<dim>::direct:
				{
		/*
		SparseDirectUMFPACK A_direct;
		A_direct.initialize(system_matrix);
		A_direct.vmult(newton_update, right_hand_side);
		*/
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
/*
        SolverControl solver_control(1000, 1e-6);
        SolverBicgstab<> solver(solver_control);

        solver.solve(system_matrix, newton_update, right_hand_side,
                    PreconditionIdentity());
        return std::pair<unsigned int, double> (solver_control.last_step(),
        		solver_control.last_value());
*/
            Assert(false,ExcNotImplemented());
            return std::pair<unsigned int, double> (0,0);
	}
	}

	// throw exception if case was not found
    Assert (false, ExcNotImplemented());
    return std::pair<unsigned int, double> (0,0);
}

template <int dim>
void ConservationLaw<dim>::compute_refinement_indicators (Vector<double> &indicator) const
{}

template <int dim>
void ConservationLaw<dim>::refine_grid (const Vector<double> &indicator)
{}

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
