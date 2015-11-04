/**
 * \file PostProcessor.cc
 * \brief Provides the function definitions for the Postprocessor class.
 */

/**
 * Constructor.
 *
 * \param[in] solution_aux_postprocessor postprocessor for derived quantities to
 *            be output
 */
template <int dim>
PostProcessor<dim>::PostProcessor(
  const ConservationLawParameters<dim> & parameters_,
  const double & end_time_,
  const bool has_exact_solution_,
  std::shared_ptr<Function<dim>> & exact_solution_function_,
  const std::string & problem_name_,
  const std::vector<std::string> & solution_component_names_,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
    solution_component_interpretations_,
  const Triangulation<dim> & triangulation_,
  const std::shared_ptr<DataPostprocessor<dim>> solution_aux_postprocessor_)
  : parameters(parameters_),
    end_time(end_time_),
    problem_name(problem_name_),
    has_exact_solution(has_exact_solution_),
    exact_solution_function(exact_solution_function_),
    solution_component_names(solution_component_names_),
    solution_component_interpretations(solution_component_interpretations_),
    is_steady_state(parameters.time_discretization ==
                    ConservationLawParameters<dim>::TemporalDiscretization::SS),
    fe(FE_Q<dim>(parameters.degree), parameters_.n_components),
    cell_quadrature(parameters.n_quadrature_points),
    current_cycle(0),
    is_last_cycle(false),
    fine_dof_handler(fine_triangulation),
    transient_output_size(0),
    transient_solution_file_number(0),
    transient_viscosity_file_number(0),
    transient_counter(0),
    transient_not_output_this_step(true),
    solution_aux_postprocessor(solution_aux_postprocessor_)
{
  // assert that a nonempty problem name was provided
  Assert(!problem_name.empty(), ExcInvalidState());

  // create map of temporal discretization to string identifier
  std::map<typename ConservationLawParameters<dim>::TemporalDiscretization,
           std::string> temp_discretization_map = {
    {ConservationLawParameters<dim>::TemporalDiscretization::ERK1, "ERK1"},
    {ConservationLawParameters<dim>::TemporalDiscretization::ERK2, "ERK2"},
    {ConservationLawParameters<dim>::TemporalDiscretization::ERK3, "ERK3"},
    {ConservationLawParameters<dim>::TemporalDiscretization::ERK4, "ERK4"},
    {ConservationLawParameters<dim>::TemporalDiscretization::SDIRK22, "SDIRK22"}};

  // determine time discretization string
  if (temp_discretization_map.find(parameters.time_discretization) ==
      temp_discretization_map.end())
  {
    ExcNotImplemented();
  }
  else
  {
    timedisc_string = temp_discretization_map[parameters.time_discretization];
  }

  // determine viscosity string
  std::string viscosity_string;
  switch (parameters.viscosity_type)
  {
    case ConservationLawParameters<dim>::ViscosityType::none:
    {
      viscosity_string = "Gal";
      break;
    }
    case ConservationLawParameters<dim>::ViscosityType::constant:
    {
      viscosity_string = "constant";
      break;
    }
    case ConservationLawParameters<dim>::ViscosityType::old_first_order:
    {
      viscosity_string = "low";
      break;
    }
    case ConservationLawParameters<dim>::ViscosityType::max_principle:
    {
      viscosity_string = "newlow";
      break;
    }
    case ConservationLawParameters<dim>::ViscosityType::entropy:
    {
      viscosity_string = "EV";
      break;
    }
    default:
    {
      ExcNotImplemented();
    }
  }

  // create filename appendage
  std::stringstream appendage_ss;
  appendage_ss << "_" << viscosity_string << "_" << timedisc_string;
  appendage_string = appendage_ss.str();

  // create filename for exact solution
  filename_exact = "solution_exact";

  // create name of output subdirectory
  std::stringstream output_dir_ss;
  output_dir_ss << "output/" << problem_name << "/";
  output_dir = output_dir_ss.str();

  // create output subdirectory
  create_directory(output_dir);

  // create fine triangulation and dof handler
  createFineTriangulationAndDoFHandler(triangulation_);

  // if outputting transient, remove older transient *.vtu files
  if (parameters.output_period > 0)
    remove_vtu_files(output_dir, "solution" + appendage_string);
}

/**
 * Destructor.
 */
template <int dim>
PostProcessor<dim>::~PostProcessor()
{
}

/**
 * \brief Outputs grid, solution, exact solution, and convergence data.
 *
 * This version has an additional parameter for a DataPostprocessor object
 * so that derived quantities can be output.
 *
 * \param[in] solution solution vector
 * \param[in] dof_handler degree of freedom handler
 * \param[in] triangulation triangulation
 */
template <int dim>
void PostProcessor<dim>::output_results(const Vector<double> & solution,
                                        const DoFHandler<dim> & dof_handler,
                                        const Triangulation<dim> & triangulation)
{
  // create output directory if it doesn't exist
  create_directory("output");

  // create output subdirectory if it doesn't exist
  create_directory(output_dir);

  // output grid
  output_grid(triangulation);

  // output solution
  std::string solution_filename = "solution" + appendage_string;
  output_solution(solution,
                  end_time,
                  dof_handler,
                  solution_filename,
                  false /* output_1d_vtu */);

  // output transient solution
  if (transient_not_output_this_step)
    output_solution_transient(solution, end_time, dof_handler, "solution", true);

  // output exact solution
  output_exact_solution(end_time);

  // output convergence data
  output_convergence_data();
}

/**
 * \brief Outputs the solution to a file.
 *
 * This version has an additional parameter for a DataPostprocessor object
 * so that derived quantities can be output.
 *
 * \param[in] values a vector of values for the quantity.
 * \param[in] time time value
 * \param[in] dof_handler degrees of freedom handler.
 * \param[in] output_string string for the output filename.
 * \param[in] output_1d_vtu option to output 1-D data in VTU format instead
 *            of GNUPLOT format
 */
template <int dim>
void PostProcessor<dim>::output_solution(const Vector<double> & solution,
                                         const double & time,
                                         const DoFHandler<dim> & dof_handler,
                                         const std::string & output_string,
                                         const bool & output_1d_vtu)
{
  // call function that outputs a vector of values to a file,
  // with the solution component names and types lists
  output_at_dof_points(solution,
                       time,
                       solution_component_names,
                       solution_component_interpretations,
                       dof_handler,
                       output_string,
                       output_1d_vtu,
                       true);
}

/**
 * \brief Outputs the solution to a file if the user specified.
 *
 * The user supplies an output period for the transient; this function
 * determines if this solution is scheduled to be output and outputs it
 * if it does.
 *
 * \param[in] solution solution vector
 * \param[in] time time value
 * \param[in] dof_handler degrees of freedom handler.
 * \param[in] output_string string for the output filename.
 * \param[in] force_output option to force output if it is not scheduled.
 */
template <int dim>
void PostProcessor<dim>::output_solution_transient(
  const Vector<double> & solution,
  const double & time,
  const DoFHandler<dim> & dof_handler,
  const std::string & output_string,
  const bool & force_output)
{
  if (parameters.output_period > 0)
  {
    // determine if this solution is scheduled to be output based on the user-
    // specified output period for the transient
    bool output_is_scheduled;
    if (transient_counter % parameters.output_period == 0)
      output_is_scheduled = true;
    else
      output_is_scheduled = false;

    // output transient solution if it is scheduled or being forced
    if (output_is_scheduled || force_output)
    {
      // make sure that transient counter is not too big for string format
      Assert(transient_counter < 10000, ExcTooManyTransientOutputFiles());

      // make sure that total transient output size has not exceeded limit
      if (transient_output_size < parameters.max_transient_output_size)
      {
        // increment transient output size estimate
        transient_output_size += 30 * solution.size();

        // create transient filename
        const std::string transient_appendage =
          "-" + Utilities::int_to_string(transient_solution_file_number, 4);

        // call function that outputs a vector of values to a file,
        // with the solution component names and types lists
        output_at_dof_points(solution,
                             time,
                             solution_component_names,
                             solution_component_interpretations,
                             dof_handler,
                             output_string + appendage_string,
                             true /* output_1d_vtu */,
                             true /* is_solution */,
                             transient_appendage);

        // signal that transient was output this step
        transient_not_output_this_step = false;

        // increment transient file number
        transient_solution_file_number++;
      }
      else
        std::cout
          << "Solution transient not output because total size limit exceeded."
          << std::endl;
    }
    else
    {
      // signal that transient was not output this step
      transient_not_output_this_step = true;
    }

    // increment transient counter
    transient_counter++;
  }
}

/**
 * \brief Outputs the viscosity to a file if the user specified.
 *
 * The user supplies an output period for the transient; this function
 * determines if this viscosity is scheduled to be output and outputs it
 * if it does.
 *
 * \param[in] solution solution vector
 * \param[in] time time value
 * \param[in] dof_handler degrees of freedom handler.
 * \param[in] force_output option to force output if it is not scheduled.
 */
template <int dim>
void PostProcessor<dim>::output_viscosity_transient(
  const std::vector<const CellMap *> & cell_maps,
  const std::vector<std::string> & names,
  const double & time,
  const DoFHandler<dim> & dof_handler,
  const bool & force_output)
{
  if (parameters.output_period > 0)
  {
    // determine if this solution is scheduled to be output based on the user-
    // specified output period for the transient
    bool output_is_scheduled;
    if (transient_counter % parameters.output_period == 0)
      output_is_scheduled = true;
    else
      output_is_scheduled = false;

    // output transient solution if it is scheduled or being forced
    if (output_is_scheduled || force_output)
    {
      // make sure that transient counter is not too big for string format
      Assert(transient_counter < 10000, ExcTooManyTransientOutputFiles());

      // make sure that total transient output size has not exceeded limit
      if (transient_output_size < parameters.max_transient_output_size)
      {
        // increment transient output size estimate
        transient_output_size += 30 * cell_maps[0]->size();

        // create transient filename
        const std::string transient_appendage =
          "-" + Utilities::int_to_string(transient_viscosity_file_number, 4);

        // call function that outputs a vector of values to a file,
        // with the solution component names and types lists
        output_cell_maps(cell_maps,
                         names,
                         "viscosity",
                         time,
                         dof_handler,
                         true /* output_1d_vtu */,
                         transient_appendage);

        // increment transient file number
        transient_viscosity_file_number++;
      }
      else
        std::cout
          << "Viscosity transient not output because total size limit exceeded."
          << std::endl;
    }
  }
}

/**
 * \brief Outputs the exact solution to a file.
 *
 * \param[in] time time value
 */
template <int dim>
void PostProcessor<dim>::output_exact_solution(const double & time)
{
  if (parameters.output_exact_solution and has_exact_solution)
  {
    // call function to evaluate and output exact solution
    output_function(*exact_solution_function,
                    time,
                    solution_component_names,
                    solution_component_interpretations,
                    filename_exact,
                    true /* is_solution */);
  }
}

/**
 * \brief Outputs to a file a vector of values, possibly with multiple
 *        components, evaluated at positions of the dofs.
 *
 * This version has an additional parameter for a DataPostprocessor object
 * so that derived quantities can be output.
 *
 * \param[in] values a vector of values for the quantity.
 * \param[in] time time value
 * \param[in] component_names list of names of each component in the
 *            vector of values
 * \param[in] component_interpretations list of types (scalar or vector)
 *            of each component in the vector of values
 * \param[in] dof_handler degrees of freedom handler.
 * \param[in] output_string string for the output filename.
 * \param[in] output_1d_vtu option to output 1-D data in VTU format instead
 *            of GNUPLOT format
 * \param[in] is_solution flag to signal that the output quantity is the
 *            solution
 * \param[in] transient_appendage string to be added onto the filename for
 *            transient files
 */
template <int dim>
void PostProcessor<dim>::output_at_dof_points(
  const Vector<double> & values,
  const double & time,
  const std::vector<std::string> & component_names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
    component_interpretations,
  const DoFHandler<dim> & dof_handler,
  const std::string & output_string,
  const bool & output_1d_vtu,
  const bool & is_solution,
  const std::string & transient_appendage)
{
  // create output directory and subdirectory if they do not exist
  create_directory("output");
  create_directory(output_dir);

  // create DataOut object for quantity
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(values,
                           component_names,
                           DataOut<dim>::type_dof_data,
                           component_interpretations);
  if (is_solution)
    if (solution_aux_postprocessor != nullptr)
      data_out.add_data_vector(values, *solution_aux_postprocessor);
  data_out.build_patches();

  // determine whether output format will be *.vtu or *.gpl
  bool output_vtu;
  if (dim > 1 || output_1d_vtu)
    output_vtu = true;
  else
    output_vtu = false;

  if (output_vtu)
  {
    // create output filestream
    std::string filename = output_string + transient_appendage + ".vtu";
    std::stringstream filename_path_ss;
    filename_path_ss << output_dir << filename;
    std::string filename_path = filename_path_ss.str();
    std::ofstream output_filestream(filename_path.c_str());
    output_filestream.precision(15);

    // write vtu file
    data_out.write_vtu(output_filestream);

    // write pvd file for associating time value to index
    if (transient_appendage != "")
    {
      times_and_solution_filenames.push_back(
        std::pair<double, std::string>(time, filename));
      std::string pvd_filename = output_dir + output_string + ".pvd";
      std::ofstream pvd_filestream(pvd_filename);
      data_out.write_pvd_record(pvd_filestream, times_and_solution_filenames);
    }
  }
  else
  {
    // create output filestream
    std::stringstream filename_ss;
    filename_ss << output_dir << output_string << ".gpl";
    std::string filename = filename_ss.str();
    std::ofstream output_filestream(filename.c_str());
    output_filestream.precision(15);

    // write output file
    data_out.write_gnuplot(output_filestream);
  }
}

/**
 * \brief Outputs convergence table to console and to a file.
 */
template <int dim>
void PostProcessor<dim>::output_convergence_data()
{
  if (has_exact_solution)
  {
    // format columns
    convergence_table.set_precision("dx", 3);
    convergence_table.set_scientific("dx", true);
    if (not is_steady_state)
    {
      convergence_table.set_precision("dt", 3);
      convergence_table.set_scientific("dt", true);
      convergence_table.set_scientific("1/dt", true);
    }
    convergence_table.set_precision("L1 error", 3);
    convergence_table.set_scientific("L1 error", true);
    convergence_table.set_precision("L2 error", 3);
    convergence_table.set_scientific("L2 error", true);

    // evaluate convergence rates
    switch (parameters.refinement_mode)
    {
      case ConservationLawParameters<dim>::RefinementMode::time:
      {
        // evaluate temporal convergence rates
        convergence_table.evaluate_convergence_rates(
          "L1 error", "1/dt", ConvergenceTable::reduction_rate_log2, 1);
        convergence_table.evaluate_convergence_rates(
          "L2 error", "1/dt", ConvergenceTable::reduction_rate_log2, 1);
        break;
      }
      case ConservationLawParameters<dim>::RefinementMode::space:
      {
        // evaluate spatial convergence rates
        convergence_table.evaluate_convergence_rates(
          "L1 error", ConvergenceTable::reduction_rate_log2);
        convergence_table.evaluate_convergence_rates(
          "L2 error", ConvergenceTable::reduction_rate_log2);
        break;
      }
      default:
      {
        ExcNotImplemented();
      }
    }

    // print convergence table to console
    std::cout << std::endl;
    convergence_table.write_text(std::cout);

    // save convergence results to file
    if (parameters.save_convergence_results)
    {
      // create output filestream for exact solution
      std::string filename =
        output_dir + "convergence" + appendage_string + ".gpl";
      std::ofstream output_filestream(filename.c_str());

      // write convergence results to file
      convergence_table.write_text(
        output_filestream, TableHandler::table_with_separate_column_description);
    }
  }
}

/**
 * \brief Outputs multiple cell maps to file.
 *
 * \param[in] cell_maps vector of pointers to cell maps
 * \param[in] names vector of names of quantities in cell maps vector
 * \param[in] filename_base base file name for output file
 * \param[in] time time value
 * \param[in] dof_handler degrees of freedom handler
 * \param[in] output_1d_vtu option to output 1-D data in VTU format instead
 *            of GNUPLOT format
 * \param[in] transient_appendage string to be added onto the filename for
 *            transient files
 */
template <int dim>
void PostProcessor<dim>::output_cell_maps(
  const std::vector<const CellMap *> & cell_maps,
  const std::vector<std::string> & names,
  const std::string & filename_base,
  const double & time,
  const DoFHandler<dim> & dof_handler,
  const bool & output_1d_vtu,
  const std::string & transient_appendage)
{
  // create output directory if it does not exist
  create_directory("output");

  // create output subdirectory if it does not exist
  create_directory(output_dir);

  // create data out object and attach DoF handler
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  // get number of cells from size of cell map
  const unsigned int n_cells = cell_maps[0]->size();

  // loop over cell maps
  const unsigned int n_cell_maps = cell_maps.size();
  std::vector<Vector<double>> cell_vectors(n_cell_maps, Vector<double>(n_cells));
  for (unsigned int i = 0; i < n_cell_maps; ++i)
  {
    // make sure the number of entries in each cell map is the same
    Assert(cell_maps[i]->size() == n_cells,
           ExcSizesInconsistent(cell_maps[i]->size(), n_cells));

    // create data vector from map
    typename CellMap::const_iterator it;
    typename CellMap::const_iterator it_end = cell_maps[i]->end();
    unsigned int j = 0;
    for (it = cell_maps[i]->begin(); it != it_end; ++it, ++j)
      cell_vectors[i][j] = it->second;

    // add data vector
    data_out.add_data_vector(
      cell_vectors[i], names[i], DataOut<dim>::type_cell_data);
  }

  // build patches
  data_out.build_patches();

  // determine whether output format will be *.vtu or *.gpl
  bool output_vtu;
  if (dim > 1 || output_1d_vtu)
    output_vtu = true;
  else
    output_vtu = false;

  // determine output file extension
  std::string filename_extension;
  if (output_vtu)
    filename_extension = ".vtu";
  else
    filename_extension = ".gpl";

  // create output filestream
  std::string filename =
    filename_base + appendage_string + transient_appendage + filename_extension;
  std::string full_filename = output_dir + filename;
  std::ofstream output_filestream(full_filename.c_str());

  // output files
  if (output_vtu)
  {
    // output *.vtu file
    data_out.write_vtu(output_filestream);

    // write pvd file for associating time value to index
    if (transient_appendage != "")
    {
      times_and_viscosity_filenames.push_back(
        std::pair<double, std::string>(time, filename));
      std::string pvd_filename =
        output_dir + filename_base + appendage_string + ".pvd";
      std::ofstream pvd_filestream(pvd_filename);
      data_out.write_pvd_record(pvd_filestream, times_and_viscosity_filenames);
    }
  }
  else
    // output *.gpl file
    data_out.write_gnuplot(output_filestream);
}

/**
 * \brief evaluate error between numerical and exact solution
 */
template <int dim>
void PostProcessor<dim>::evaluate_error(const Vector<double> & solution,
                                        const DoFHandler<dim> & dof_handler,
                                        const Triangulation<dim> & triangulation)
{
  if (has_exact_solution)
  {
    // assert that this function is only being called when an exact solution is
    // available
    Assert(has_exact_solution, ExcInvalidState());

    // set time for exact solution function
    exact_solution_function->set_time(end_time);

    // number of cells
    unsigned int n_cells = triangulation.n_active_cells();

    // error per cell
    Vector<double> difference_per_cell(n_cells);

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
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    double domain_volume = 0.0;
    for (; cell != endc; ++cell)
      domain_volume += cell->measure();
    double avg_cell_volume = domain_volume / n_cells;

    // add error values to convergence table
    convergence_table.add_value("cycle", current_cycle);
    convergence_table.add_value("cells", n_cells);
    convergence_table.add_value("dofs", dof_handler.n_dofs());
    convergence_table.add_value("dx", avg_cell_volume);
    if (not is_steady_state)
    {
      convergence_table.add_value("dt", dt_nominal);
      convergence_table.add_value("1/dt", 1.0 / dt_nominal);
    }
    convergence_table.add_value("L1 error", L1_error);
    convergence_table.add_value("L2 error", L2_error);
  }
}

/**
 * \brief Outputs the grid.
 *
 * \param[in] triangulation the triangulation to output
 */
template <int dim>
void PostProcessor<dim>::output_grid(
  const Triangulation<dim> & triangulation) const
{
  if (parameters.output_mesh)
  {
    if (dim > 1)
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
    else
    {
      std::cout
        << "User specified to output the mesh, but 1-d meshes cannot be output"
        << std::endl;
    }
  }
}

/**
 * \brief Updates the time step size to be put in convergence table
 *
 * \param[in] dt_nominal nominal time step size for this cycle
 */
template <int dim>
void PostProcessor<dim>::update_dt(const double & dt)
{
  dt_nominal = dt;
}

/**
 * \brief Checks if a directory exists and creates it if it does not
 *
 * \param[in] directory path to directory to create if it does not exist
 */
template <int dim>
void PostProcessor<dim>::create_directory(const std::string & directory) const
{
  // convert to char
  char * directory_char = (char *)directory.c_str();

  // use stat to determine if directory exists
  struct stat mystat;
  bool directory_exists = false;
  if (stat(directory_char, &mystat) == 0)
    if (S_ISDIR(mystat.st_mode))
      directory_exists = true;

  // create directory if it doesn't exist
  int make_status;
  if (!directory_exists)
    make_status = system(("mkdir " + directory).c_str());
  else
    make_status = 0;
  Assert(make_status == 0, ExcInternalError());
}

/**
 * \brief Sets the current cycle and flags it if it is the last.
 *
 * \param[in] cycle current refinement cycle
 */
template <int dim>
void PostProcessor<dim>::set_cycle(const unsigned int & cycle)
{
  current_cycle = cycle;
  is_last_cycle = (cycle == parameters.n_refinement_cycles - 1);
}

/**
 * \brief Creates fine triangulation and dof handler for outputting functions.
 *
 * \param[in] triangulation triangulation to copy for refining
 */
template <int dim>
void PostProcessor<dim>::createFineTriangulationAndDoFHandler(
  const Triangulation<dim> & triangulation)
{
  // create fine mesh on which to interpolate functions
  fine_triangulation.copy_triangulation(triangulation);
  const unsigned int final_refinement_level =
    parameters.initial_refinement_level + parameters.n_refinement_cycles - 1;
  const int n_refinements =
    parameters.exact_solution_refinement_level - final_refinement_level;
  if (n_refinements > 0)
    fine_triangulation.refine_global(n_refinements);
}

/**
 * \brief Computes a function on a fine triangulation and outputs to file.
 *
 * \param[in] function the function to be evaluated
 * \param[in] time time value
 * \param[in] component_names list of names of each component in the
 *            vector of values
 * \param[in] component_interpretations list of types (scalar or vector)
 *            of each component in the vector of values
 * \param[in] filename the output filename
 * \param[in] is_solution flag to signal that the output quantity is the
 *            solution
 */
template <int dim>
void PostProcessor<dim>::output_function(
  Function<dim> & function,
  const double & time,
  const std::vector<std::string> & component_names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
    component_interpretations,
  const std::string & filename,
  const bool & is_solution)
{
  // create FESystem object for function
  const unsigned int n_components_function = component_names.size();
  FESystem<dim> fe_function(FE_Q<dim>(parameters.degree), n_components_function);

  // create dof handler for function and distribute dofs
  DoFHandler<dim> dof_handler(fine_triangulation);
  dof_handler.distribute_dofs(fe_function);

  // compute function values
  Vector<double> function_values(dof_handler.n_dofs());
  function.set_time(time);
  VectorTools::interpolate(dof_handler, function, function_values);

  // output function values to file
  output_at_dof_points(function_values,
                       time,
                       component_names,
                       component_interpretations,
                       dof_handler,
                       filename,
                       false /* output_1d_vtu */,
                       is_solution);
}

/**
 * \brief Removes all filename_base*.vtu files in a directory.
 *
 * \param[in] directory directory for which to delete .vtu files
 * \param[in] filename_base filename base; any file matching filename*.vtu
 *            is removed
 */
template <int dim>
void PostProcessor<dim>::remove_vtu_files(const std::string & directory,
                                          const std::string & filename_base) const
{
  // create a regular expression for *.vtu files
  std::regex vtu_regex(filename_base + ".*\\.vtu");

  // use dirent.h to read a directory's contents
  DIR * dir;
  struct dirent * ent;

  // check to make sure directory can be opened
  AssertThrow((dir = opendir(directory.c_str())) != NULL,
              ExcDirectoryCannotBeOpened(directory));

  // loop through files in directory
  while ((ent = readdir(dir)) != NULL)
  {
    // get filename
    std::string filename = ent->d_name;

    // determine if file is a *.vtu file
    if (regex_match(filename, vtu_regex))
    {
      // create full file path
      std::string filepath = directory + filename;

      // remove the *.vtu file
      std::remove(filepath.c_str());
    }
  }
  // close the directory
  closedir(dir);
}
