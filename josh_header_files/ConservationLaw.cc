/*
template <int dim>
void ConservationLaw<dim>::print_test_message_from_base_class()
{
  std::cout << "successfully printed message" << std::endl;
}
*/

template <int dim>
ConservationLaw<dim>::ConservationLaw(ParameterHandler &prm,//const std::string &input_filename,
                                      const int &n_comp):
//   conservation_law_parameters(n_comp),
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

   assemble_system();

   // output initial solution
   output_results();

   // begin time stepping loop
   double time = 0;
   double next_time_step_output = conservation_law_parameters.output_period;
   while (time < conservation_law_parameters.final_time)
   {
	   std::cout << "time: " << time << std::endl
			   << "   Number of active cells:       "
			   << triangulation.n_active_cells()
			   << std::endl
			   << "   Number of degrees of freedom: "
			   << dof_handler.n_dofs()
			   << std::endl
			   << std::endl;

	   std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl
			   << "   _____________________________________" << std::endl;

	   unsigned int nonlin_iter = 0;
	   current_solution = predictor;
	   // initialize convergence flag
	   bool not_converged_yet = true;
	   // begin Newton loop for current time step
	   while (not_converged_yet)
	   {
		   system_matrix = 0;
		   right_hand_side = 0;
		   assemble_system ();

		   const double res_norm = right_hand_side.l2_norm();
		   if (std::fabs(res_norm) < 1e-10)
		   {
			   // nonlinear solution converged
			   std::printf("   %-16.3e (converged)\n\n", res_norm);
			   break;
		   }
		   else
		   {
			   // nonlinear solution has not yet converged; take another Newton step
			   Vector<double> newton_update = 0;

			   std::pair<unsigned int, double> convergence
			   = solve (newton_update);

			   current_solution += newton_update;

			   std::printf("   %-16.3e %04d        %-5.2e\n",
					   res_norm, convergence.first, convergence.second);
		   }

		   ++nonlin_iter;
		   AssertThrow (nonlin_iter <= 10,
				   ExcMessage ("No convergence in nonlinear solver"));
	   }
	   time += conservation_law_parameters.time_step;

	   if (conservation_law_parameters.output_period < 0)
		   output_results ();
	   else if (time >= next_time_step_output)
	   {
		   output_results ();
		   next_time_step_output += conservation_law_parameters.output_period;
	   }

	   predictor = current_solution;
	   predictor.sadd (2.0, -1.0, old_solution);

	   old_solution = current_solution;

	   /*
	   if (parameters.do_refine == true)
	   {
		   Vector<double> refinement_indicators (triangulation.n_active_cells());
		   compute_refinement_indicators(refinement_indicators);

		   refine_grid(refinement_indicators);
		   setup_system();

		   newton_update.reinit (dof_handler.n_dofs());
	   }
	   */
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
   std::cout << "Setup system" << std::endl;
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
std::pair<unsigned int, double> ConservationLaw<dim>::solve (Vector<double> &solution)
{
   std::pair<unsigned int, double> solutions;
   std::cout << "Solve" << std::endl;
   return solutions;
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
	//typename EulerEquations<dim>::Postprocessor
	//postprocessor (parameters.schlieren_plot);

	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);

	data_out.add_data_vector (current_solution,
			component_names,
			DataOut<dim>::type_dof_data,
			component_interpretations);

	data_out.add_data_vector (current_solution, "current_solution");//postprocessor);

	data_out.build_patches ();

	static unsigned int output_file_number = 0;
	std::string filename = "solution-" +
			Utilities::int_to_string (output_file_number, 3) +
			".vtk";
	std::ofstream output (filename.c_str());
	data_out.write_vtk (output);

	++output_file_number;
}
