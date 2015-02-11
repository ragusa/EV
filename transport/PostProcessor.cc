/** \brief Constructor.
 */
template<int dim>
PostProcessor<dim>::PostProcessor(
   const bool          &output_mesh,
   const bool          &output_exact_solution,
   const bool          &save_convergence_results,
   const bool          &has_exact_solution,
   FunctionParser<dim> &exact_solution_function,
   const double        &time,
   const double        &dt_nominal,
   const bool          &is_steady_state,
   const unsigned int  &refinement_option,
   const unsigned int  &final_refinement_level,
   const FESystem<dim> &fe,
   const unsigned int  &degree,
   const unsigned int  &scheme_option,
   const unsigned int  &problem_ID,
   const QGauss<dim>   &cell_quadrature) :
   output_mesh(output_mesh),
   output_exact_solution(output_exact_solution),
   save_convergence_results(save_convergence_results),
   has_exact_solution(has_exact_solution),
   exact_solution_function(&exact_solution_function),
   time(time),
   dt_nominal(dt_nominal),
   is_steady_state(is_steady_state),
   refinement_option(refinement_option),
   final_refinement_level(final_refinement_level),
   fe(&fe),
   degree(degree),
   scheme_option(scheme_option),
   problem_ID(problem_ID),
   cell_quadrature(cell_quadrature)
{
   // determine viscosity string
   switch (scheme_option) {
      case 0: {
         viscosity_string = "_galerkin";
         break;
      } case 1: {
         viscosity_string = "_low_order";
         break;
      } case 2: {
         viscosity_string = "_high_order";
         break;
      } case 3: {
         viscosity_string = "_FCT";
         break;
      } default: {
         Assert(false, ExcNotImplemented());
         break;
      }
   }
}

/** \brief Destructor.
 */
template<int dim>
PostProcessor<dim>::~PostProcessor()
{
}

/** \brief Output grid, solution, and viscosity to output file and print
 *         convergence table.
 */
template<int dim>
void PostProcessor<dim>::output_results(const Vector<double>     &solution,
                                        const DoFHandler<dim>    &dof_handler,
                                        const Triangulation<dim> &triangulation)
{
   // output grid
   //------------
   if ((output_mesh) and (dim > 1))
      output_grid(triangulation);

   // output solution
   output_solution(solution,
                   dof_handler,
                   "solution",
                   true);

   // write exact solution output file
   //---------------------------------
   if (output_exact_solution and has_exact_solution)
   {
      // create fine mesh on which to interpolate exact solution function
      Triangulation<dim> fine_triangulation;
      fine_triangulation.copy_triangulation(triangulation);
      // define "fine" triangulation to be of refinement level 10, so refine
      // if refinement level is below this
      if (final_refinement_level < 10) fine_triangulation.refine_global(10-final_refinement_level);
      // create dof handler for fine mesh
      DoFHandler<dim> fine_dof_handler(fine_triangulation);
      fine_dof_handler.distribute_dofs(*fe);
      // interpolate exact solution
      exact_solution_function->set_time(time);
      Vector<double> exact_solution(fine_dof_handler.n_dofs());
      VectorTools::interpolate(fine_dof_handler,
                               *exact_solution_function,
                               exact_solution);

      // output exact solution to file
      output_solution(exact_solution,
                      fine_dof_handler,
                      "exact_solution",
                      false);
   }


/*
   // output min and max bounds for DMP
   if (parameters.output_DMP_bounds and parameters.scheme_option == 3)
   {
      output_solution(min_values,dof_handler,"min_values",false);
      output_solution(max_values,dof_handler,"max_values",false);
   }
*/

   // output convergence table
   //------------------------
   if (has_exact_solution) {
      // format columns
      convergence_table.set_precision("dx", 3);
      convergence_table.set_scientific("dx", true);
      if (is_steady_state) {
         convergence_table.set_precision("dt", 3);
         convergence_table.set_scientific("dt", true);
         convergence_table.set_scientific("1/dt", true);
      }
      convergence_table.set_precision("L1 error", 3);
      convergence_table.set_scientific("L1 error", true);
      convergence_table.set_precision("L2 error", 3);
      convergence_table.set_scientific("L2 error", true);
      if (refinement_option == 2) {
         // evaluate temporal convergence rates
         convergence_table.evaluate_convergence_rates("L1 error", "1/dt", ConvergenceTable::reduction_rate_log2, 1);
         convergence_table.evaluate_convergence_rates("L2 error", "1/dt", ConvergenceTable::reduction_rate_log2, 1);
      } else {
         // evaluate spatial convergence rates
         convergence_table.evaluate_convergence_rates("L1 error", ConvergenceTable::reduction_rate_log2);
         convergence_table.evaluate_convergence_rates("L2 error", ConvergenceTable::reduction_rate_log2);
      }
      // print convergence table to console
      std::cout << std::endl;
      convergence_table.write_text(std::cout);

      // save convergence results to file
      if (save_convergence_results) {
         // create filename
         std::stringstream filename_ss;
         filename_ss << "output/convergence" << viscosity_string << "_" << problem_ID << ".gpl";
         std::string filename = filename_ss.str();
         char *filename_char = (char*)filename.c_str();
         // create output filestream for exact solution
         std::ofstream output_filestream(filename_char);
         // write convergence results to file
         convergence_table.write_text(output_filestream,TableHandler::table_with_separate_column_description);
      }
   }
}

/** \brief Outputs a solution to a file.
 *  \param [in] solution a vector of solution values.
 *  \param [in] dof_handler degrees of freedom handler.
 *  \param [in] output_string string for the output filename.
 *  \param [in] append_viscosity the option to include a string for viscosity type in output filename
 */
template<int dim>
void PostProcessor<dim>::output_solution(const Vector<double>  &solution,
                                         const DoFHandler<dim> &dof_handler,
                                         const std::string     &output_string,
                                         const bool            &append_viscosity) const
{
   // create DataOut object for solution
   DataOut<dim> data_out;
   data_out.attach_dof_handler(dof_handler);
   data_out.add_data_vector(solution, "flux");
   data_out.build_patches(degree + 1);

   // create string for viscosity type if it is to be included in output filename
   std::string append_string;
   if (append_viscosity)
      append_string = viscosity_string;
   else
      append_string = "";

   // create output filename for exact solution
   std::string filename_extension;
   if (dim == 1) filename_extension = ".gpl";
   else          filename_extension = ".vtk";

   std::stringstream filename_ss;
   filename_ss << "output/" << output_string << append_string
      << "_" << problem_ID << filename_extension;
   std::string filename = filename_ss.str();
   char *filename_char = (char*)filename.c_str();

   // create output filestream for exact solution
   std::ofstream output_filestream(filename_char);
   // write file
   if (dim == 1) data_out.write_gnuplot(output_filestream);
   else          data_out.write_vtk    (output_filestream);
}

/** \brief Outputs viscosities to file.
 *  \param [in] low_order_viscosity low-order viscosity in each cell.
 *  \param [in] entropy_viscosity entropy viscosity in each cell.
 *  \param [in] high_order_viscosity high-order viscosity in each cell.
 *  \param [in] dof_handler degrees of freedom handler.
 */
template<int dim>
void PostProcessor<dim>::output_viscosity(const Vector<double>  &low_order_viscosity,
                                          const Vector<double>  &entropy_viscosity,
                                          const Vector<double>  &high_order_viscosity,
                                          const DoFHandler<dim> &dof_handler) const
{
   // write viscosity output file
   //----------------------------
   if (scheme_option != 0) {
      DataOut<dim> visc_out;
      visc_out.attach_dof_handler(dof_handler);
      // add viscosity data vector(s)
      switch (scheme_option)
      {
         case 1: {
            visc_out.add_data_vector(low_order_viscosity,"Low_Order_Viscosity",DataOut<dim>::type_cell_data);
            break;
         } case 2: {
            visc_out.add_data_vector(low_order_viscosity, "Low_Order_Viscosity", DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(entropy_viscosity,   "Entropy_Viscosity",   DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(high_order_viscosity,"High_Order_Viscosity",DataOut<dim>::type_cell_data);
            break;
         } case 3: {
            visc_out.add_data_vector(low_order_viscosity, "Low_Order_Viscosity", DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(entropy_viscosity,   "Entropy_Viscosity",   DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(high_order_viscosity,"High_Order_Viscosity",DataOut<dim>::type_cell_data);
            break;
         } default: {
            Assert(false, ExcNotImplemented());
            break;
         }
      }

      // determine output file extension
      std::string filename_extension;
      if (dim ==  1)
         filename_extension = ".gpl";
      else
         filename_extension = ".vtk";

      // create output filename
      std::stringstream viscosity_filename_ss;
      viscosity_filename_ss << "output/viscosity" << filename_extension;
      std::string viscosity_filename = viscosity_filename_ss.str();
      char *viscosity_filename_char = (char*)viscosity_filename.c_str();

      // create output filestream
      std::ofstream viscosity_outstream(viscosity_filename_char);
      // build patches and write to file
      visc_out.build_patches(degree + 1);
      if (dim == 1)
         visc_out.write_gnuplot(viscosity_outstream);
      else
         visc_out.write_vtk(viscosity_outstream);
   }
}

/** \brief evaluate error between numerical and exact solution
 */
template<int dim>
void PostProcessor<dim>::evaluate_error(const Vector<double>     &solution,
                                        const DoFHandler<dim>    &dof_handler,
                                        const Triangulation<dim> &triangulation,
                                        const unsigned int       &cycle)
{
   // assert that this function is only being called when an exact solution is available
   Assert(has_exact_solution,ExcInvalidState());

   // set time for exact solution function
   exact_solution_function->set_time(time);

   // number of cells
   unsigned int n_cells = triangulation.n_active_cells();

   // error per cell
   Vector<double> difference_per_cell (n_cells);

   // compute L1 error
   VectorTools::integrate_difference(MappingQ<dim>(1),
                                     dof_handler,
                                     solution,
                                     *exact_solution_function,
                                     difference_per_cell,
                                     cell_quadrature,
                                     VectorTools::L1_norm);
   const double L1_error = difference_per_cell.l1_norm();

   // compute L2 error
   VectorTools::integrate_difference(MappingQ<dim>(1),
                                     dof_handler,
                                     solution,
                                     *exact_solution_function,
                                     difference_per_cell,
                                     cell_quadrature,
                                     VectorTools::L2_norm);
   const double L2_error = difference_per_cell.l2_norm();

   // compute average cell volume
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   double domain_volume = 0.0;
   for (; cell != endc; ++cell)
      domain_volume += cell->measure();
   double avg_cell_volume = domain_volume / n_cells;

   // add error values to convergence table
   convergence_table.add_value("cycle", cycle);
   convergence_table.add_value("cells", n_cells);
   convergence_table.add_value("dofs", dof_handler.n_dofs());
   convergence_table.add_value("dx", avg_cell_volume);
   if (not is_steady_state) {
      convergence_table.add_value("dt", dt_nominal);
      convergence_table.add_value("1/dt", 1.0/dt_nominal);
   }
   convergence_table.add_value("L1 error", L1_error);
   convergence_table.add_value("L2 error", L2_error);
}

/** \brief output the grid of the given cycle
 */
template<int dim>
void PostProcessor<dim>::output_grid(const Triangulation<dim> &triangulation) const
{
   // create output filestream
   std::string filename = "grid.eps";
   std::ofstream output(filename.c_str());

   // write grid to eps
   GridOut grid_out;
   grid_out.write_eps(triangulation, output);
}
