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
   const typename RefinementHandler<dim>::RefinementMode &refinement_mode,
   const unsigned int  &final_refinement_level,
   const FESystem<dim> &fe,
   const std::string   &output_dir,
   const std::string   &appendage_string,
   const std::string   &filename_exact,
   const QGauss<dim>   &cell_quadrature) :

   output_mesh(output_mesh),
   output_exact_solution(output_exact_solution),
   save_convergence_results(save_convergence_results),
   has_exact_solution(has_exact_solution),
   exact_solution_function(&exact_solution_function),
   time(time),
   dt_nominal(dt_nominal),
   is_steady_state(is_steady_state),
   refinement_mode(refinement_mode),
   final_refinement_level(final_refinement_level),
   fe(&fe),
   output_dir(output_dir),
   appendage_string(appendage_string),
   filename_exact(filename_exact),
   cell_quadrature(cell_quadrature)
{
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
   // create output directory if it doesn't exist
   create_directory("output");

   // create output subdirectory if it doesn't exist
   create_directory(output_dir);

   // output grid
   if (output_mesh) {
      if (dim > 1) output_grid(triangulation);
      else std::cout << "User specified to output the mesh, but 1-d meshes cannot be outputted" << std::endl;
   }

   // output solution
   std::string solution_filename = "solution" + appendage_string;
   output_solution(solution,
                   dof_handler,
                   solution_filename);

   // write exact solution output file
   if (output_exact_solution and has_exact_solution)
   {
      // create fine mesh on which to interpolate exact solution function
      Triangulation<dim> fine_triangulation;
      fine_triangulation.copy_triangulation(triangulation);
      // define "fine" triangulation to be of refinement level 7, so refine
      // if refinement level is below this
      unsigned int fine_level = 7;
      if (final_refinement_level < fine_level) fine_triangulation.refine_global(
         fine_level-final_refinement_level);
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
                      filename_exact);
   }

   // output convergence table
   if (has_exact_solution) {
      // format columns
      convergence_table.set_precision("dx", 3);
      convergence_table.set_scientific("dx", true);
      if (not is_steady_state) {
         convergence_table.set_precision("dt", 3);
         convergence_table.set_scientific("dt", true);
         convergence_table.set_scientific("1/dt", true);
      }
      convergence_table.set_precision("L1 error", 3);
      convergence_table.set_scientific("L1 error", true);
      convergence_table.set_precision("L2 error", 3);
      convergence_table.set_scientific("L2 error", true);
      switch (refinement_mode) {
         case RefinementHandler<dim>::time : {
            // evaluate temporal convergence rates
            convergence_table.evaluate_convergence_rates("L1 error", "1/dt", ConvergenceTable::reduction_rate_log2, 1);
            convergence_table.evaluate_convergence_rates("L2 error", "1/dt", ConvergenceTable::reduction_rate_log2, 1);
            break;
         } case RefinementHandler<dim>::space : {
            // evaluate spatial convergence rates
            convergence_table.evaluate_convergence_rates("L1 error", ConvergenceTable::reduction_rate_log2);
            convergence_table.evaluate_convergence_rates("L2 error", ConvergenceTable::reduction_rate_log2);
            break;
         } default : {
            ExcNotImplemented();
         }
      }
      // print convergence table to console
      std::cout << std::endl;
      convergence_table.write_text(std::cout);

      // save convergence results to file
      if (save_convergence_results) {
         // create output filestream for exact solution
         std::string filename = output_dir + "convergence" + appendage_string + ".gpl";
         std::ofstream output_filestream(filename.c_str());
         // write convergence results to file
         convergence_table.write_text(output_filestream,TableHandler::table_with_separate_column_description);
      }
   }
}

/** \brief Outputs a solution to a file.
 *  \param [in] solution a vector of solution values.
 *  \param [in] dof_handler degrees of freedom handler.
 *  \param [in] output_string string for the output filename.
 */
template<int dim>
void PostProcessor<dim>::output_solution(const Vector<double>  &solution,
                                         const DoFHandler<dim> &dof_handler,
                                         const std::string     &output_string) const
{
   // create output directory if it doesn't exist
   create_directory("output");

   // create output subdirectory if it doesn't exist
   create_directory(output_dir);

   // create DataOut object for solution
   DataOut<dim> data_out;
   data_out.attach_dof_handler(dof_handler);
   data_out.add_data_vector(solution, "flux");
   data_out.build_patches();

   // create output filename for solution
   std::string filename_extension;
   if (dim == 1) filename_extension = ".gpl";
   else          filename_extension = ".vtk";

   std::stringstream filename_ss;
   filename_ss << output_dir << output_string << filename_extension;
   std::string filename = filename_ss.str();

   // create output filestream for exact solution
   std::ofstream output_filestream(filename.c_str());
   output_filestream.precision(15);
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
   // create output directory if it doesn't exist
   create_directory("output");

   // create output subdirectory if it doesn't exist
   create_directory(output_dir);

   // add viscosities to data out object
   DataOut<dim> visc_out;
   visc_out.attach_dof_handler(dof_handler);
   visc_out.add_data_vector(low_order_viscosity, "Low_Order_Viscosity", DataOut<dim>::type_cell_data);
   visc_out.add_data_vector(entropy_viscosity,   "Entropy_Viscosity",   DataOut<dim>::type_cell_data);
   visc_out.add_data_vector(high_order_viscosity,"High_Order_Viscosity",DataOut<dim>::type_cell_data);

   // determine output file extension
   std::string filename_extension;
   if (dim ==  1) filename_extension = ".gpl";
   else           filename_extension = ".vtk";

   // create output filestream
   std::string viscosity_filename = output_dir + "viscosity" + appendage_string
      + filename_extension;
   std::ofstream viscosity_outstream(viscosity_filename.c_str());

   // build patches and write to file
   visc_out.build_patches();
   if (dim == 1) visc_out.write_gnuplot(viscosity_outstream);
   else          visc_out.write_vtk    (viscosity_outstream);
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
   // create output directory if it doesn't exist
   create_directory("output");

   // create output subdirectory if it doesn't exist
   create_directory(output_dir);

   // create output filestream
   std::string filename = output_dir + "grid.eps";
   std::ofstream output(filename.c_str());

   // write grid to eps
   GridOut grid_out;
   grid_out.write_eps(triangulation, output);
}

/** \brief Update the time step size to be put in convergence table
 */
template<int dim>
void PostProcessor<dim>::update_dt(const double &dt)
{
   dt_nominal = dt;
}

/** \brief Check if a directory exists and create it if it doesn't.
 */
template<int dim>
void PostProcessor<dim>::create_directory(const std::string &directory) const
{
   // convert to char
   char *directory_char = (char*)directory.c_str();

   // use stat to determine if directory exists
   struct stat mystat;
   bool directory_exists = false;
   if (stat(directory_char,&mystat) == 0)
      if (S_ISDIR(mystat.st_mode))
         directory_exists = true;

   // create directory if it doesn't exist
   int make_status = 0;
   if (!directory_exists)
      make_status = system(("mkdir "+directory).c_str());
   Assert(make_status == 0, ExcInternalError());
}
