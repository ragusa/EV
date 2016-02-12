/** \brief constructor for TransportProblem class
 */
template <int dim>
TransportProblem<dim>::TransportProblem(
  const TransportParameters<dim> & parameters)
  : mpi_communicator(MPI_COMM_WORLD),
    parameters(parameters),
    triangulation(mpi_communicator),
    dof_handler(triangulation),
    degree(parameters.degree),
    fe(degree),
    dofs_per_cell(fe.dofs_per_cell),
    faces_per_cell(GeometryInfo<dim>::faces_per_cell),
    cell_quadrature(parameters.n_quadrature_points),
    face_quadrature(parameters.n_quadrature_points),
    n_q_points_cell(cell_quadrature.size()),
    n_q_points_face(face_quadrature.size()),
    has_exact_solution(false),
    domain_volume(0.0),
    pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    computing_timer(
      mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times)
{
}

/** \brief destructor for TransportProblem class
 */
template <int dim>
TransportProblem<dim>::~TransportProblem()
{
  dof_handler.clear();
}

/** \brief initialize system
 */
template <int dim>
void TransportProblem<dim>::initialize_system()
{
  // timer
  TimerOutput::Scope t_initialize(computing_timer, "initialize");

  // process problem ID
  process_problem_ID();

  // decide that transport direction is the unit x vector
  transport_direction = 0;
  transport_direction[0] = 1.0;

  // determine function variables based on dimension;
  // this is to be used in initialization of function parser objects
  std::string variables;
  if (dim == 1)
    variables = "x";
  else if (dim == 2)
    variables = "x,y";
  else if (dim == 3)
    variables = "x,y,z";
  else
    Assert(false, ExcInvalidState());

  // create function parser constants
  function_parser_constants["pi"] = numbers::PI;
  function_parser_constants["x_min"] = x_min;
  function_parser_constants["x_mid"] = x_min + 0.5 * (x_max - x_min);
  function_parser_constants["x_max"] = x_max;

  // initialize exact solution function
  if (has_exact_solution)
    exact_solution_function.initialize(
      variables, exact_solution_string, function_parser_constants, false);

  // initialize source function
  source_function.initialize(
    variables, source_string, function_parser_constants, false);

  // initialize cross section function
  cross_section_function.initialize(
    variables, cross_section_string, function_parser_constants, false);

  // initialize Dirichlet boundary value function
  incoming_function.initialize(
    variables, incoming_string, function_parser_constants, false);

  // create grid for initial refinement level
  GridGenerator::hyper_cube(triangulation, x_min, x_max);
  domain_volume = std::pow((x_max - x_min), dim);
  triangulation.refine_global(parameters.initial_refinement_level);
  n_cells = triangulation.n_active_cells();
}

/** \brief process problem ID
 */
template <int dim>
void TransportProblem<dim>::process_problem_ID()
{
  switch (parameters.problem_id)
  {
    case 1:
    { // pure absorber
      Assert(dim < 3, ExcNotImplemented());

      x_min = 0.0;
      x_max = 1.0;

      incoming_string = "1";
      function_parser_constants["incoming"] = 1.0;

      cross_section_string = "1.0";
      function_parser_constants["sigma"] = 1.0;

      source_string = "0";
      function_parser_constants["source"] = 0.0;

      has_exact_solution = true;
      exact_solution_string =
        "source/sigma + (incoming - source/sigma)*exp(-sigma*(x-x_min))";

      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }
}

/** \brief set up the problem before assembly of the linear system
 */
template <int dim>
void TransportProblem<dim>::setup_system()
{
  // timer
  TimerOutput::Scope t_setup(computing_timer, "setup");

  // distribute dofs
  dof_handler.distribute_dofs(fe);
  n_dofs = dof_handler.n_dofs();

  // get index sets of locally owned and relevant dofs
  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  // reinitialize solution vector, system matrix, and rhs
  solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  // set boundary indicators to distinguish incoming boundary
  set_boundary_indicators();

  // clear constraint matrix and make hanging node constraints
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(
    dof_handler, 1, incoming_function, constraints);
  constraints.close();

  // create sparsity pattern for system matrix and mass matrices
  CompressedSimpleSparsityPattern csp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false);
  SparsityTools::distribute_sparsity_pattern(
    csp,
    dof_handler.n_locally_owned_dofs_per_processor(),
    mpi_communicator,
    locally_relevant_dofs);

  // reinitialize system matrix and mass matrices
  system_matrix.reinit(
    locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator);
}

/** \brief Assemble the inviscid system matrix. The inviscid steady-state matrix
 *         is independent of the solution and independent of time and thus needs
 *         to be called only once per level of mesh refinement.
 */
template <int dim>
void TransportProblem<dim>::assemble_system()
{
  // timer
  TimerOutput::Scope t_assemble_matrix(computing_timer, "assemble");

  // FE values, for assembly terms
  FEValues<dim> fe_values(fe,
                          cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // total cross section and source values at each quadrature point on cell
  std::vector<double> total_cross_section_values(n_q_points_cell);
  std::vector<double> source_values(n_q_points_cell);

  // cell iterator
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  // loop over cells
  for (; cell != endc; ++cell)
  {
    if (cell->is_locally_owned())
    {
      // initialize local matrix and rhs to zero
      cell_matrix = 0;
      cell_rhs = 0;

      // reinitialize FE values
      fe_values.reinit(cell);

      // get quadrature points on cell
      std::vector<Point<dim>> points(n_q_points_cell);
      points = fe_values.get_quadrature_points();

      // get cross section values for all quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
      {
        total_cross_section_values[q] = cross_section_function.value(points[q]);
        source_values[q] = source_function.value(points[q]);
      }

      // compute cell contributions to global system
      // ------------------------------------------------------------------
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int q = 0; q < n_q_points_cell; ++q)
      {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          cell_rhs(i) +=
            fe_values.shape_value(i, q) * source_values[q] * fe_values.JxW(q);
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            cell_matrix(i, j) +=
              (
                // divergence term
                fe_values.shape_value(i, q) * transport_direction *
                  fe_values.shape_grad(j, q) +
                // total interaction term
                fe_values.shape_value(i, q) * total_cross_section_values[q] *
                  fe_values.shape_value(j, q)) *
              fe_values.JxW(q);
          } // end j
        }   // end i
      }     // end q

      // aggregate local matrix and rhs to global matrix and rhs
      constraints.distribute_local_to_global(
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    }
  } // end cell

  // compress matrix
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

/** \brief Sets the boundary indicators for each boundary face.
 *
 *         The Dirichlet BC is applied only to the incoming boundary, so the
 * transport
 *         direction is compared against the normal vector of the face.
 */
template <int dim>
void TransportProblem<dim>::set_boundary_indicators()
{
  // cell iterator
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  // reset boundary indicators to zero
  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
    if (cell->is_locally_owned())
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary())
          cell->face(face)->set_boundary_indicator(0);

  // FE face values
  FEFaceValues<dim> fe_face_values(fe, face_quadrature, update_normal_vectors);
  // loop over cells
  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
  {
    if (cell->is_locally_owned())
    {
      // loop over faces of cell
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
           ++face)
      {
        // if face is at boundary
        if (cell->face(face)->at_boundary())
        {
          // reinitialize FE face values
          fe_face_values.reinit(cell, face);
          // determine if the transport flux is incoming through this face;
          //  it isn't necessary to loop over all face quadrature points because
          //  the transport direction and normal vector are the same at each
          //  quadrature point; therefore, quadrature point 0 is arbitrarily
          //  chosen
          double small = -1.0e-12;
          if (fe_face_values.normal_vector(0) * transport_direction < small)
          {
            // mark boundary as incoming flux boundary: indicator 1
            cell->face(face)->set_boundary_indicator(1);
          }
        }
      }
    }
  }
}

/** \brief run the problem
 */
template <int dim>
void TransportProblem<dim>::run()
{
  // loop over refinement cycles
  for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; ++cycle)
  {
    if (cycle == 0)
      // initialize problem
      initialize_system();
    else
      // refine mesh
      refine_grid();

    // setup system
    setup_system();

    // print information
    pcout << std::endl << "Cycle " << cycle << ':' << std::endl;
    pcout << "   Number of active cells:       " << n_cells << std::endl;
    pcout << "   Number of degrees of freedom: " << n_dofs << std::endl;

    // assemble system
    assemble_system();

    // solve linear system
    solve();

    // output solution
    if (cycle == parameters.n_refinement_cycles - 1)
      output_solution();

    // print timer summary and reset
    computing_timer.print_summary();
    computing_timer.reset();
  }

  // prevent timer from continuing to print output
  computing_timer.disable_output();
}

/** \brief Solves the steady-state system.
 */
template <int dim>
void TransportProblem<dim>::solve()
{
  // timer
  TimerOutput::Scope t_solve(computing_timer, "solve");

  LA::MPI::Vector completely_distributed_solution(locally_owned_dofs,
                                                  mpi_communicator);

  // linear solver tolerance
  const double linear_tol = 1e-10 * system_rhs.l2_norm() + 1e-10;
  SolverControl solver_control(dof_handler.n_dofs(), linear_tol);

  // preconditioner
  LA::MPI::PreconditionJacobi preconditioner;
  preconditioner.initialize(system_matrix);

  // solve
  PETScWrappers::SolverBiCG solver(solver_control, mpi_communicator);
  solver.solve(
    system_matrix, completely_distributed_solution, system_rhs, preconditioner);

  pcout << "Solved in " << solver_control.last_step() << " iterations."
        << std::endl;

  constraints.distribute(completely_distributed_solution);

  solution = completely_distributed_solution;
}

template <int dim>
void TransportProblem<dim>::output_solution()
{
  // timer
  TimerOutput::Scope t_output(computing_timer, "output");

  if (parameters.output_solution)
  {
    // write to vtk file, so assume not 1-D
    Assert(dim > 1, ExcNotImplemented());

    // create DataOut object for solution
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "flux");

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches();

    // create output filename for solution
    std::string filename =
      ("solution-" +
       Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4));

    // create output filestream and write file
    std::ofstream output((filename + ".vtu").c_str());
    data_out.write_vtu(output);

    // write master record
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      // create list of filenames for each subdomain
      std::vector<std::string> filenames;
      for (unsigned int i = 0;
           i < Utilities::MPI::n_mpi_processes(mpi_communicator);
           ++i)
        filenames.push_back("solution-" + Utilities::int_to_string(i, 4) +
                            ".vtu");
      // create output filestream for master record and write to file
      std::ofstream master_output("solution.pvtu");
      data_out.write_pvtu_record(master_output, filenames);
    }
  }
}

/** \brief Refine the grid.
 */
template <int dim>
void TransportProblem<dim>::refine_grid()
{
  Assert(parameters.use_adaptive_refinement == false, ExcNotImplemented());

  // refine uniformly
  triangulation.refine_global(1);

  // update number of cells
  n_cells = triangulation.n_active_cells();
}
