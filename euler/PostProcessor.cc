/**
 * \file PostProcessor.cc
 * \brief Provides the function definitions for the Postprocessor class.
 */

/**
 * Constructor.
 */
template <int dim>
PostProcessor<dim>::PostProcessor(
  const ConservationLawParameters<dim> & parameters_,
  const bool has_exact_solution_,
  std::shared_ptr<Function<dim>> & exact_solution_function_,
  const std::string & problem_name_,
  const std::vector<std::string> & solution_component_names_,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
    solution_component_interpretations_,
  const Triangulation<dim> & triangulation_)
  : parameters(parameters_),
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
    fine_dof_handler(fine_triangulation)
{
  // assert that a nonempty problem name was provided
  Assert(!problem_name.empty(), ExcInvalidState());

  // create map of temporal discretization to string identifier
  std::map<typename ConservationLawParameters<dim>::TemporalDiscretization,
           std::string> temp_discretization_map = {
    /*
        {ConservationLawParameters<dim>::TemporalDiscretization::SS,"SS"},
        {ConservationLawParameters<dim>::TemporalDiscretization::FE,"FE"},
        {ConservationLawParameters<dim>::TemporalDiscretization::CN,"CN"},
        {ConservationLawParameters<dim>::TemporalDiscretization::BE,"BE"},
        {ConservationLawParameters<dim>::TemporalDiscretization::SSP2,"SSPRK22"},
        {ConservationLawParameters<dim>::TemporalDiscretization::SSP3,"SSPRK33"}
    */
    {ConservationLawParameters<dim>::TemporalDiscretization::ERK1, "ERK1"},
    {ConservationLawParameters<dim>::TemporalDiscretization::ERK2, "ERK2"},
    {ConservationLawParameters<dim>::TemporalDiscretization::ERK3, "ERK3"},
    {ConservationLawParameters<dim>::TemporalDiscretization::ERK4, "ERK4"},
    {ConservationLawParameters<dim>::TemporalDiscretization::SDIRK22, "SDIRK22"}};

  // determine time discretization string
  std::string timedisc_string;
  if (temp_discretization_map.find(parameters.time_discretization) ==
      temp_discretization_map.end())
  {
    ExcNotImplemented();
  }
  else
  {
    timedisc_string = temp_discretization_map[parameters.time_discretization];
  }

  /*
    switch (parameters.time_discretization_option)
    {
      case ConservationLawParameters<dim>::TemporalDiscretization::SS :
      {
        timedisc_string = "SS";
        break;
      }
      case ConservationLawParameters<dim>::TemporalDiscretization::FE :
      {
        timedisc_string = "FE";
        break;
      }
      case ConservationLawParameters<dim>::TemporalDiscretization::CN :
      {
        timedisc_string = "CN";
        break;
      }
      case ConservationLawParameters<dim>::TemporalDiscretization::BE :
      {
        timedisc_string = "BE";
        break;
      }
      case ConservationLawParameters<dim>::TemporalDiscretization::SSP2 :
      {
        timedisc_string = "SSPRK22";
        break;
      }
      case ConservationLawParameters<dim>::TemporalDiscretization::SSP3 :
      {
        timedisc_string = "SSPRK33";
        break;
      }
      default :
      {
        ExcNotImplemented();
      }
    }
  */

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

  // create fine triangulation and dof handler
  createFineTriangulationAndDoFHandler(triangulation_);
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
  output_solution(solution, dof_handler, solution_filename);

  // output exact solution
  output_exact_solution();

  // output convergence data
  output_convergence_data();
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
 * \param[in] data_postprocessor postprocessor for derived quantities to
 *            be output
 */
template <int dim>
void PostProcessor<dim>::output_results(const Vector<double> & solution,
                                        const DoFHandler<dim> & dof_handler,
                                        const Triangulation<dim> & triangulation,
                                        const DataPostprocessor<dim> & data_postprocessor)
{
  // create output directory if it doesn't exist
  create_directory("output");

  // create output subdirectory if it doesn't exist
  create_directory(output_dir);

  // output grid
  output_grid(triangulation);

  // output solution
  std::string solution_filename = "solution" + appendage_string;
  output_solution(solution, dof_handler, solution_filename, data_postprocessor);

  // output exact solution
  output_exact_solution();

  // output convergence data
  output_convergence_data();
}

/**
 * \brief Outputs the solution to a file.
 *
 * \param[in] values a vector of values for the quantity.
 * \param[in] dof_handler degrees of freedom handler.
 * \param[in] output_string string for the output filename.
 */
template <int dim>
void PostProcessor<dim>::output_solution(const Vector<double> & solution,
                                         const DoFHandler<dim> & dof_handler,
                                         const std::string & output_string) const
{
  // call function that outputs a vector of values to a file,
  // with the solution component names and types lists
  output_at_dof_points(solution,
                       solution_component_names,
                       solution_component_interpretations,
                       dof_handler,
                       output_string);
}

/**
 * \brief Outputs the solution to a file.
 *
 * This version has an additional parameter for a DataPostprocessor object
 * so that derived quantities can be output.
 *
 * \param[in] values a vector of values for the quantity.
 * \param[in] dof_handler degrees of freedom handler.
 * \param[in] output_string string for the output filename.
 * \param[in] data_postprocessor postprocessor for derived quantities to
 *            be output
 */
template <int dim>
void PostProcessor<dim>::output_solution(const Vector<double> & solution,
                                         const DoFHandler<dim> & dof_handler,
                                         const std::string & output_string,
                                         const DataPostprocessor<dim> & data_postprocessor) const
{
  // call function that outputs a vector of values to a file,
  // with the solution component names and types lists
  output_at_dof_points(solution,
                       solution_component_names,
                       solution_component_interpretations,
                       dof_handler,
                       output_string,
                       data_postprocessor);
}

/**
 * \brief Outputs the exact solution to a file.
 */
template <int dim>
void PostProcessor<dim>::output_exact_solution()
{
  if (parameters.output_exact_solution and has_exact_solution)
  {
    // call function to evaluate and output exact solution
    output_function(*exact_solution_function,
                    solution_component_names,
                    solution_component_interpretations,
                    filename_exact,
                    parameters.end_time);
  }
}

/**
 * \brief Outputs to a file a vector of values, possibly with multiple
 *        components, evaluated at positions of the dofs.
 *
 * \param[in] values a vector of values for the quantity.
 * \param[in] component_names list of names of each component in the
 *            vector of values
 * \param[in] component_interpretations list of types (scalar or vector)
 *            of each component in the vector of values
 * \param[in] dof_handler degrees of freedom handler.
 * \param[in] output_string string for the output filename.
 */
template <int dim>
void PostProcessor<dim>::output_at_dof_points(
  const Vector<double> & values,
  const std::vector<std::string> & component_names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
    component_interpretations,
  const DoFHandler<dim> & dof_handler,
  const std::string & output_string) const
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
    data_out.build_patches();

    // write GNUplot file for 1-D and .vtk file otherwise
    if (dim == 1)
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
    else
    {
      // create output filestream
      std::stringstream filename_ss;
      filename_ss << output_dir << output_string << ".vtk";
      std::string filename = filename_ss.str();
      std::ofstream output_filestream(filename.c_str());
      output_filestream.precision(15);

      // write output file
      data_out.write_vtk(output_filestream);
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
 * \param[in] component_names list of names of each component in the
 *            vector of values
 * \param[in] component_interpretations list of types (scalar or vector)
 *            of each component in the vector of values
 * \param[in] dof_handler degrees of freedom handler.
 * \param[in] output_string string for the output filename.
 * \param[in] data_postprocessor postprocessor for derived quantities to
 *            be output
 */
template <int dim>
void PostProcessor<dim>::output_at_dof_points(
  const Vector<double> & values,
  const std::vector<std::string> & component_names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
    component_interpretations,
  const DoFHandler<dim> & dof_handler,
  const std::string & output_string,
  const DataPostprocessor<dim> & data_postprocessor) const
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
    data_out.add_data_vector(values, data_postprocessor);
    data_out.build_patches();

    // write GNUplot file for 1-D and .vtk file otherwise
    if (dim == 1)
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
    else
    {
      // create output filestream
      std::stringstream filename_ss;
      filename_ss << output_dir << output_string << ".vtk";
      std::string filename = filename_ss.str();
      std::ofstream output_filestream(filename.c_str());
      output_filestream.precision(15);

      // write output file
      data_out.write_vtk(output_filestream);
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

/** \brief Outputs viscosities to file.
 *  \param [in] low_order_viscosity low-order viscosity in each cell.
 *  \param [in] entropy_viscosity entropy viscosity in each cell.
 *  \param [in] high_order_viscosity high-order viscosity in each cell.
 *  \param [in] dof_handler degrees of freedom handler.
 */
template <int dim>
void PostProcessor<dim>::output_viscosity(
  const Vector<double> & low_order_viscosity,
  const Vector<double> & entropy_viscosity,
  const Vector<double> & high_order_viscosity,
  const DoFHandler<dim> & dof_handler) const
{
  // create output directory if it doesn't exist
  create_directory("output");

  // create output subdirectory if it doesn't exist
  create_directory(output_dir);

  // add viscosities to data out object
  DataOut<dim> visc_out;
  visc_out.attach_dof_handler(dof_handler);
  visc_out.add_data_vector(
    low_order_viscosity, "Low_Order_Viscosity", DataOut<dim>::type_cell_data);
  visc_out.add_data_vector(
    entropy_viscosity, "Entropy_Viscosity", DataOut<dim>::type_cell_data);
  visc_out.add_data_vector(
    high_order_viscosity, "High_Order_Viscosity", DataOut<dim>::type_cell_data);

  // determine output file extension
  std::string filename_extension;
  if (dim == 1)
    filename_extension = ".gpl";
  else
    filename_extension = ".vtk";

  // create output filestream
  std::string viscosity_filename =
    output_dir + "viscosity" + appendage_string + filename_extension;
  std::ofstream viscosity_outstream(viscosity_filename.c_str());

  // build patches and write to file
  visc_out.build_patches();
  if (dim == 1)
    visc_out.write_gnuplot(viscosity_outstream);
  else
    visc_out.write_vtk(viscosity_outstream);
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
    exact_solution_function->set_time(parameters.end_time);

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

/** \brief output the grid of the given cycle
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

/** \brief Update the time step size to be put in convergence table
 */
template <int dim>
void PostProcessor<dim>::update_dt(const double & dt)
{
  dt_nominal = dt;
}

/** \brief Check if a directory exists and create it if it doesn't.
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
  int make_status = 0;
  if (!directory_exists)
    make_status = system(("mkdir " + directory).c_str());
  Assert(make_status == 0, ExcInternalError());
}

/**
 * \brief Sets the current cycle and flags it if it is the last.
 */
template <int dim>
void PostProcessor<dim>::set_cycle(const unsigned int & cycle)
{
  current_cycle = cycle;
  is_last_cycle = (cycle == parameters.n_refinement_cycles - 1);
}

/**
 * \brief Creates fine triangulation and dof handler for outputting functions.
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
 * This version is for time-independent functions.
 *
 * \param[in] function the function to be evaluated
 * \param[in] component_names list of names of each component in the
 *            vector of values
 * \param[in] component_interpretations list of types (scalar or vector)
 *            of each component in the vector of values
 * \param[in] filename the output filename
 */
template <int dim>
void PostProcessor<dim>::output_function(
  const Function<dim> & function,
  const std::vector<std::string> & component_names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
    component_interpretations,
  const std::string & filename) const
{
  // create FESystem object for function
  const unsigned int n_components_function = component_names.size();
  FESystem<dim> fe_function(FE_Q<dim>(parameters.degree), n_components_function);

  // create dof handler for function and distribute dofs
  DoFHandler<dim> dof_handler(fine_triangulation);
  dof_handler.distribute_dofs(fe_function);

  // compute function values
  Vector<double> function_values(dof_handler.n_dofs());
  VectorTools::interpolate(dof_handler, function, function_values);

  // output function values to file
  output_at_dof_points(function_values,
                       component_names,
                       component_interpretations,
                       dof_handler,
                       filename);
}

/**
 * \brief Computes a function on a fine triangulation and outputs to file.
 *
 * This version is for time-dependent functions.
 *
 * \param[in] function the function to be evaluated
 * \param[in] component_names list of names of each component in the
 *            vector of values
 * \param[in] component_interpretations list of types (scalar or vector)
 *            of each component in the vector of values
 * \param[in] filename the output filename
 * \param[in] t the time at which to evaluate the function
 */
template <int dim>
void PostProcessor<dim>::output_function(
  Function<dim> & function,
  const std::vector<std::string> & component_names,
  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &
    component_interpretations,
  const std::string & filename,
  const double & t) const
{
  // set time of function
  function.set_time(t);

  // call time-independent version of function
  output_function(function, component_names, component_interpretations, filename);
}
