/**
 * \file ProblemParameters.cc
 * \brief Provides the function definitions for the ProblemParameters class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] problem_name_            name of problem
 * \param[in] n_components_            number of solution components
 * \param[in] specified_steady_state_  flag that problem is to be run in
 *            steady-state
 */
template <int dim>
ProblemParameters<dim>::ProblemParameters(const std::string & problem_name_,
                                          const unsigned int & n_components_,
                                          const bool & specified_steady_state_)
  : problem_name(problem_name_),
    n_components(n_components_),
    specified_steady_state(specified_steady_state_),
    initial_conditions_function(n_components_)
{
}

/**
 * \brief Declares base problem parameters.
 */
template <int dim>
void ProblemParameters<dim>::declare_base_parameters()
{
  // validity
  parameter_handler.enter_subsection("validity");
  {
    parameter_handler.declare_entry("valid in 1d",
                                    "false",
                                    Patterns::Bool(),
                                    "Flag signalling problem is valid in 1-D");
    parameter_handler.declare_entry("valid in 2d",
                                    "false",
                                    Patterns::Bool(),
                                    "Flag signalling problem is valid in 2-D");
    parameter_handler.declare_entry("valid in 3d",
                                    "false",
                                    Patterns::Bool(),
                                    "Flag signalling problem is valid in 3-D");
    parameter_handler.declare_entry(
      "is transient problem",
      "true",
      Patterns::Bool(),
      "Flag signalling that problem is a transient problem only");
  }
  parameter_handler.leave_subsection();

  // domain
  parameter_handler.enter_subsection("domain");
  {
    parameter_handler.declare_entry("domain shape",
                                    "hyper_cube",
                                    Patterns::Anything(),
                                    "Descriptor for shape of domain");
    parameter_handler.declare_entry(
      "x start", "0.0", Patterns::Double(), "Start of domain in x-direction");
    parameter_handler.declare_entry(
      "y start", "0.0", Patterns::Double(), "Start of domain in y-direction");
    parameter_handler.declare_entry(
      "z start", "0.0", Patterns::Double(), "Start of domain in z-direction");
    parameter_handler.declare_entry(
      "x width", "1.0", Patterns::Double(), "Width of domain in x-direction");
    parameter_handler.declare_entry(
      "y width", "1.0", Patterns::Double(), "Width of domain in y-direction");
    parameter_handler.declare_entry(
      "z width", "1.0", Patterns::Double(), "Width of domain in z-direction");
  }
  parameter_handler.leave_subsection();

  // boundary conditions
  parameter_handler.enter_subsection("boundary conditions");
  {
    parameter_handler.declare_entry("boundary conditions type",
                                    "none",
                                    Patterns::Anything(),
                                    "Type of boundary conditions");
    parameter_handler.declare_entry("use exact solution as dirichlet bc",
                                    "false",
                                    Patterns::Bool(),
                                    "Flag signalling that exact solution is to"
                                    "be used as Dirichlet BC");
    parameter_handler.declare_entry(
      "boundary id scheme",
      "all",
      Patterns::Selection("all|incoming|left"),
      "Option for how to assign boundary ID to faces");
  }
  parameter_handler.leave_subsection();

  // initial conditions
  parameter_handler.enter_subsection("initial conditions");
  {
    parameter_handler.declare_entry("use exact solution as initial conditions",
                                    "false",
                                    Patterns::Bool(),
                                    "Flag signalling that exact solution is to"
                                    "be used as initial conditions");
  }
  parameter_handler.leave_subsection();

  // exact solution
  parameter_handler.enter_subsection("exact solution");
  {
    parameter_handler.declare_entry("has exact solution",
                                    "false",
                                    Patterns::Bool(),
                                    "Flag signalling that exact solution exists");
    parameter_handler.declare_entry("exact solution type",
                                    "function",
                                    Patterns::Anything(),
                                    "Type of exact solution");
  }
  parameter_handler.leave_subsection();

  // default end time
  parameter_handler.enter_subsection("default end time");
  {
    parameter_handler.declare_entry(
      "has default end time",
      "false",
      Patterns::Bool(),
      "Flag signalling that a default end time exists");
    parameter_handler.declare_entry(
      "default end time", "1.0", Patterns::Double(), "Default end time");
  }
  parameter_handler.leave_subsection();
}

/**
 * \brief Gets and processes base and derived problem parameters.
 *
 * \param[in] parameters_file  parameters file name
 * \param[in] triangulation    triangulation
 * \param[in] fe               finite element system
 * \param[in] face_quadrature  face quadrature
 */
template <int dim>
void ProblemParameters<dim>::get_and_process_parameters(
  const std::string & parameters_file,
  Triangulation<dim> & triangulation,
  const FESystem<dim> & fe,
  const QGauss<dim - 1> & face_quadrature)
{
  // declare base and derived parameters
  declare_base_parameters();
  declare_derived_parameters();

  // assert that parameters input file exists
  struct stat buffer;
  const bool file_exists = stat(parameters_file.c_str(), &buffer) == 0;
  Assert(file_exists, ExcFileDoesNotExist(parameters_file));

  // read parameters input file
  parameter_handler.read_input(parameters_file);

  // get base and derived parameters
  get_base_parameters();
  get_derived_parameters();

  // generate initial mesh, compute domain volume, and set boundary IDs
  generate_mesh_and_compute_volume(triangulation);
  set_boundary_ids(triangulation, fe, face_quadrature);

  // process derived and base parameters
  process_base_parameters();
  process_derived_parameters(triangulation, fe, face_quadrature);
  process_shared_base_parameters(triangulation, fe, face_quadrature);
}

/**
 * \brief Gets base problem parameters from parameter handler.
 */
template <int dim>
void ProblemParameters<dim>::get_base_parameters()
{
  // validity
  parameter_handler.enter_subsection("validity");
  {
    valid_in_1d = parameter_handler.get_bool("valid in 1d");
    valid_in_2d = parameter_handler.get_bool("valid in 2d");
    valid_in_3d = parameter_handler.get_bool("valid in 3d");
    is_transient_problem = parameter_handler.get_bool("is transient problem");
  }
  parameter_handler.leave_subsection();

  // domain
  parameter_handler.enter_subsection("domain");
  {
    domain_shape = parameter_handler.get("domain shape");
    x_start = parameter_handler.get_double("x start");
    y_start = parameter_handler.get_double("y start");
    z_start = parameter_handler.get_double("z start");
    x_width = parameter_handler.get_double("x width");
    y_width = parameter_handler.get_double("y width");
    z_width = parameter_handler.get_double("z width");
  }
  parameter_handler.leave_subsection();

  // initial conditions
  parameter_handler.enter_subsection("initial conditions");
  {
    use_exact_solution_as_initial_conditions =
      parameter_handler.get_bool("use exact solution as initial conditions");
  }
  parameter_handler.leave_subsection();

  // boundary conditions
  parameter_handler.enter_subsection("boundary conditions");
  {
    boundary_conditions_type = parameter_handler.get("boundary conditions type");
    use_exact_solution_as_dirichlet_bc =
      parameter_handler.get_bool("use exact solution as dirichlet bc");
    boundary_id_scheme = parameter_handler.get("boundary id scheme");
  }
  parameter_handler.leave_subsection();

  // exact solution
  parameter_handler.enter_subsection("exact solution");
  {
    has_exact_solution = parameter_handler.get_bool("has exact solution");
    exact_solution_type = parameter_handler.get("exact solution type");
  }
  parameter_handler.leave_subsection();

  // default end time
  parameter_handler.enter_subsection("default end time");
  {
    has_default_end_time = parameter_handler.get_bool("has default end time");
    default_end_time = parameter_handler.get_double("default end time");
  }
  parameter_handler.leave_subsection();
}

/**
 * \brief Processes the base parameters to provide necessary public data.
 */
template <int dim>
void ProblemParameters<dim>::process_base_parameters()
{
  // assert number of dimensions is valid
  if (!valid_in_1d)
  {
    Assert(dim != 1, ExcImpossibleInDim(dim));
  }
  if (!valid_in_2d)
  {
    Assert(dim != 2, ExcImpossibleInDim(dim));
  }
  if (!valid_in_3d)
  {
    Assert(dim != 3, ExcImpossibleInDim(dim));
  }

  // assert that problem is not transient if steady-state is specified
  Assert(!(specified_steady_state && is_transient_problem),
         ExcNotASteadyStateProblem());

  // constants for function parsers
  constants["pi"] = numbers::PI;
  constants["x_min"] = x_start;
  constants["x_max"] = x_start + x_width;
  constants["x_mid"] = x_start + 0.5 * x_width;
  constants["y_min"] = y_start;
  constants["y_max"] = y_start + y_width;
  constants["z_min"] = z_start;
  constants["z_max"] = z_start + z_width;
}

/**
 * \brief Processes the shared base parameters to provide necessary public data.
 *
 * \param[in] triangulation    triangulation
 * \param[in] fe               finite element system
 * \param[in] face_quadrature  face quadrature
 */
template <int dim>
void ProblemParameters<dim>::process_shared_base_parameters(
  Triangulation<dim> &,
  const FESystem<dim> & fe,
  const QGauss<dim - 1> & face_quadrature)
{
  // exact solution
  if (has_exact_solution)
  {
    // if exact solution was not created in derived class
    if (!(exact_solution_function))
    {
      if (exact_solution_type == "parsed")
      {
        // create and initialize function parser for exact solution
        auto exact_solution_function_derived =
          std::make_shared<FunctionParser<dim>>(n_components);
        exact_solution_function_derived->initialize(
          FunctionParser<dim>::default_variable_names() + ",t",
          exact_solution_strings,
          constants,
          true);

        // create base class function pointer
        exact_solution_function = exact_solution_function_derived;
      }
      else
      {
        AssertThrow(false, ExcNotImplemented());
      }
    }
  }

  // if boundary conditions were not created in derived class
  if (!(boundary_conditions))
  {
    // create boundary conditions
    if (boundary_conditions_type == "dirichlet")
    {
      auto derived_boundary_conditions =
        std::make_shared<DirichletBoundaryConditions<dim>>(fe, face_quadrature);
      boundary_conditions = derived_boundary_conditions;

      // initialize Dirichlet boundary functions if needed
      if (use_exact_solution_as_dirichlet_bc)
      {
        dirichlet_function = exact_solution_function;
      }
      else
      {
        auto dirichlet_function_derived =
          std::make_shared<FunctionParser<dim>>(n_components);
        dirichlet_function_derived->initialize(
          FunctionParser<dim>::default_variable_names() + ",t",
          dirichlet_function_strings,
          constants,
          true);

        // point base class pointer to derived class object
        dirichlet_function = dirichlet_function_derived;
      }
    }
    else
    {
      AssertThrow(false, ExcNotImplemented());
    }
  }

  // if chose to use exact solution as initial conditions
  if (use_exact_solution_as_initial_conditions)
    initial_conditions_strings = exact_solution_strings;

  // initialize initial conditions function
  initial_conditions_function.initialize(
    FunctionParser<dim>::default_variable_names() + ",t",
    initial_conditions_strings,
    constants,
    true);
}

/**
 * \brief Generates the mesh and computes the domain volume.
 *
 * \param[in] triangulation  triangulation
 */
template <int dim>
void ProblemParameters<dim>::generate_mesh_and_compute_volume(
  Triangulation<dim> & triangulation)
{
  if (domain_shape == "hyper_cube")
  {
    // compute domain volume
    domain_volume = std::pow(x_width, dim);

    // generate mesh
    GridGenerator::hyper_cube(triangulation, x_start, x_start + x_width);
  }
  else if (domain_shape == "hyper_box")
  {
    // diagonally opposite begin and end points defining box
    Point<dim> point_begin;
    Point<dim> point_end;

    // store constants for function parsers
    const double x_end = x_start + x_width;
    point_begin[0] = x_start;
    point_end[0] = x_end;
    constants["x_min"] = x_start;
    constants["x_max"] = x_end;
    domain_volume = x_width; // initialize domain volume

    if (dim > 1)
    {
      const double y_end = y_start + y_width;
      point_begin[1] = y_start;
      point_end[1] = y_end;
      constants["y_min"] = y_start;
      constants["y_max"] = y_end;
      domain_volume *= y_width; // update volume
    }
    else // dim == 3
    {
      const double z_end = z_start + z_width;
      point_begin[2] = z_start;
      point_end[2] = z_end;
      constants["z_min"] = z_start;
      constants["z_max"] = z_end;
      domain_volume *= z_width; // update volume
    }

    // generate mesh
    GridGenerator::hyper_rectangle(triangulation, point_begin, point_end);
  }
  else if (domain_shape == "2d_dam_break")
  {
    AssertThrow(dim == 2, ExcImpossibleInDim(dim));

    // create vertices
    std::vector<Point<2>> vertices = {
      Point<2>(0, 0), Point<2>(1, 0), Point<2>(2, 0), Point<2>(4, 0),
      Point<2>(5, 0), Point<2>(6, 0), Point<2>(0, 1), Point<2>(1, 1),
      Point<2>(2, 1), Point<2>(4, 1), Point<2>(5, 1), Point<2>(6, 1),
      Point<2>(0, 2), Point<2>(1, 2), Point<2>(2, 2), Point<2>(3, 2),
      Point<2>(4, 2), Point<2>(5, 2), Point<2>(6, 2), Point<2>(0, 3),
      Point<2>(1, 3), Point<2>(2, 3), Point<2>(3, 3), Point<2>(4, 3),
      Point<2>(5, 3), Point<2>(6, 3), Point<2>(0, 4), Point<2>(1, 4),
      Point<2>(2, 4), Point<2>(3, 4), Point<2>(4, 4), Point<2>(5, 4),
      Point<2>(6, 4), Point<2>(0, 5), Point<2>(1, 5), Point<2>(2, 5),
      Point<2>(4, 5), Point<2>(5, 5), Point<2>(6, 5), Point<2>(0, 6),
      Point<2>(1, 6), Point<2>(2, 6), Point<2>(4, 6), Point<2>(5, 6),
      Point<2>(6, 6)};

    // cells
    const unsigned int n_cells_coarse = 28;
    const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell] = {
      {0, 1, 6, 7},     {1, 2, 7, 8},     {3, 4, 9, 10},    {4, 5, 10, 11},
      {6, 7, 12, 13},   {7, 8, 13, 14},   {9, 10, 16, 17},  {10, 11, 17, 18},
      {12, 13, 19, 20}, {13, 14, 20, 21}, {14, 15, 21, 22}, {15, 16, 22, 23},
      {16, 17, 23, 24}, {17, 18, 24, 25}, {19, 20, 26, 27}, {20, 21, 27, 28},
      {21, 22, 28, 29}, {22, 23, 29, 30}, {23, 24, 30, 31}, {24, 25, 31, 32},
      {26, 27, 33, 34}, {27, 28, 34, 35}, {30, 31, 36, 37}, {31, 32, 37, 38},
      {33, 34, 39, 40}, {34, 35, 40, 41}, {36, 37, 42, 43}, {37, 38, 43, 44}};

    std::vector<CellData<2>> cell_data(n_cells_coarse, CellData<2>());
    for (unsigned int i = 0; i < n_cells_coarse; ++i)
    {
      for (unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_cell; ++j)
        cell_data[i].vertices[j] = cell_vertices[i][j];
      cell_data[i].material_id = 0;
    }

    // domain volume
    domain_volume = 28.0;

    // create triangulation
    triangulation.create_triangulation(vertices, cell_data, SubCellData());
  }
  else
  {
    Assert(false, ExcNotImplemented());
  }

  // assert that domain volume be nonzero
  Assert(std::fabs(domain_volume) > 1.0e-15, ExcInvalidState());
}

/**
 * \brief Sets the boundary IDs.
 *
 * \param[in] triangulation  triangulation
 * \param[in] fe               finite element system
 * \param[in] face_quadrature  face quadrature
 */
template <int dim>
void ProblemParameters<dim>::set_boundary_ids(
  Triangulation<dim> & triangulation,
  const FESystem<dim> & fe,
  const QGauss<dim - 1> & face_quadrature)
{
  // compute faces per cell
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

  // FE face values for getting quadrature points on faces
  FEFaceValues<dim> fe_face_values(fe, face_quadrature, update_quadrature_points);

  // loop over cells
  typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
  typename Triangulation<dim>::cell_iterator endc = triangulation.end();
  for (; cell != endc; ++cell)
    // loop over faces
    for (unsigned int face = 0; face < faces_per_cell; ++face)
      // if face is a boundary face
      if (cell->face(face)->at_boundary())
      {
        // get quadrature points on face
        fe_face_values.reinit(cell, face);
        std::vector<Point<dim>> points(fe_face_values.get_quadrature_points());

        // set boundary ID for face
        if (boundary_id_scheme == "all")
        {
          cell->face(face)->set_boundary_id(0);
        }
        else if (boundary_id_scheme == "left")
        {
          // left boundary should have ID of 0
          if (points[0][0] < x_start + 1.0e-10)
            cell->face(face)->set_boundary_id(0);
          else
            cell->face(face)->set_boundary_id(1);
        }
        else
        {
          // "incoming" boundary ID scheme is the only scheme implemented in
          // a derived class, so if the selected scheme does not match any scheme
          // in this branch, it should mean that the scheme is "incoming"
          Assert(boundary_id_scheme == "incoming", ExcNotImplemented());
        }
      }
}
