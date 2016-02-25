/**
 * \file ProblemParameters.cc
 * \brief Provides the function definitions for the ProblemParameters class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ProblemParameters<dim>::ProblemParameters()
{
}

/**
 * \brief Declares base problem parameters.
 */
template <int dim>
void ProblemParameters<dim>::declare_base_parameters()
{
  // dimension
  parameter_handler.enter_subsection("dimension");
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
 * \param[in] filename       parameters file name
 * \param[in] triangulation  triangulation
 */
template <int dim>
void ProblemParameters<dim>::get_and_process_parameters(
  const std::string & filename,
  Triangulation<dim> & triangulation)
{
  // declare base and derived parameters
  declare_base_parameters();
  declare_derived_parameters();

  // assert that parameters input file exists
  struct stat buffer;
  const bool file_exists = stat(filename.c_str(), &buffer) == 0;
  Assert(file_exists, ExcFileDoesNotExist(filename));

  // read parameters input file
  parameter_handler.read_input(filename);

  // get base and derived parameters
  get_base_parameters();
  get_derived_parameters();

  // process derived and base parameters
  process_derived_parameters();
  process_base_parameters(triangulation);
}

/**
 * \brief Gets base problem parameters from parameter handler.
 */
template <int dim>
void ProblemParameters<dim>::get_base_parameters()
{
  // dimension
  parameter_handler.enter_subsection("dimension");
  {
    valid_in_1d = parameter_handler.get_bool("valid in 1d");
    valid_in_2d = parameter_handler.get_bool("valid in 2d");
    valid_in_3d = parameter_handler.get_bool("valid in 3d");
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
 *
 * \param[in] triangulation  triangulation
 */
template <int dim>
void ProblemParameters<dim>::process_base_parameters(
  Triangulation<dim> & triangulation)
{
  // assert number of dimensions is valid
  if (!parameters.valid_in_1d)
    Assert(dim != 1, ExcImpossibleInDim(dim));
  if (!parameters.valid_in_2d)
    Assert(dim != 2, ExcImpossibleInDim(dim));
  if (!parameters.valid_in_3d)
    Assert(dim != 3, ExcImpossibleInDim(dim));

  // constants for function parsers
  constants["pi"] = numbers::PI;

  // exact solution
  if (has_exact_solution)
  {
    // if exact solution was not created in derived class
    if (!(exact_solution))
    {
      if (exact_solution_type == "function")
      {
        // create and initialize function parser for exact solution
        std::shared_ptr<FunctionParser<dim>> exact_solution_function_derived =
          std::make_shared<FunctionParser<dim>>(1);
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
        Assert(false, ExcNotImplemented());
      }
    }
  }

  // create boundary conditions
  if (boundary_conditions_type == "dirichlet")
  {
    auto derived_boundary_conditions =
      std::make_shared<DirichletBoundaryConditions<dim>>(this->fe,
                                                         this->face_quadrature);
    boundary_conditions = derived_boundary_conditions;

  // initialize Dirichlet boundary functions if needed
  if (use_exact_solution_as_dirichlet_bc)
  {
      n_dirichlet_boundaries = 1;
    dirichlet_function.resize(1);
    dirichlet_function[0] = exact_solution_function;
  }
  else
  {
    dirichlet_function.resize(n_dirichlet_boundaries);
    for (unsigned int boundary = 0; boundary < n_dirichlet_boundaries; ++boundary)
    {
      dirichlet_function[boundary] = std::make_shared<FunctionParser<dim>>(n_components);
      dirichlet_function[boundary]->initialize(
        FunctionParser<dim>::default_variable_names() + ",t",
        dirichlet_function_strings[boundary],
        constants,
        true);
    }
  }
  }
  else
  {
    Assert(false, ExcNotImplemented());
  }

  // initialize initial conditions function
  initial_conditions_function.initialize(
    FunctionParser<dim>::default_variable_names(),
    initial_conditions_strings, constants, false);

  // default end time
  has_default_end_time = parameters.has_default_end_time;
  default_end_time = parameters.default_end_time;
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
  if (parameters.domain_shape == "hyper_cube")
  {
    // compute domain volume
    const double x_start = parameters.x_start;
    const double x_width = parameters.x_width;
    domain_volume = std::pow(x_width, dim);

    // generate mesh
    GridGenerator::hyper_cube(triangulation, x_start, x_start + x_width);
  }
  else if (parameters.domain_shape == "hyper_box")
  {
    // store constants for function parsers
    const double x_width = parameters.x_width;
    constants["x_min"] = parameters.x_start;
    constants["x_max"] = parameters.x_start + x_width;
    domain_volume = x_width; // initialize domain volume

    if (dim > 1)
    {
      const double y_width = parameters.y_width;
      constants["y_min"] = parameters.y_start;
      constants["y_max"] = parameters.y_start + y_width;
      domain_volume *= y_width; // update volume
    }
    else // dim == 3
    {
      const double z_width = parameters.z_width;
      constants["z_min"] = parameters.z_start;
      constants["z_max"] = parameters.z_start + z_width;
      domain_volume *= z_width; // update volume
    }
  }
  else
  {
    Assert(false, ExcNotImplemented());
  }
}

/**
 * \brief Sets the boundary IDs.
 *
 * \param[in] triangulation  triangulation
 */
template <int dim>
void ProblemParameters<dim>::set_boundary_ids(Triangulation<dim> & triangulation)
{
  // compute faces per cell
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

  // loop over cells
  typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
  typename Triangulation<dim>::cell_iterator endc = triangulation.end();
  for (; cell != endc; ++cell)
    // loop over faces
    for (unsigned int face = 0; face < faces_per_cell; ++face)
      // if face is a boundary face
      if (cell->face(face)->at_boundary())
      {
        // set boundary ID for face
        if (boundary_id_scheme == "all")
          cell->face(face)->set_boundary_id(0);
        else
        {
          // "incoming" boundary ID scheme is the only scheme implemented in
          // a derived class, so if the selected scheme does not match any scheme
          // in this branch, it should mean that the scheme is "incoming"
          Assert(boundary_id_scheme == "incoming", ExcNotImplemented());
        }
      }
}
