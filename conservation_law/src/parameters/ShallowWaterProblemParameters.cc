/**
 * \file ShallowWaterProblemParameters.cc
 * \brief Provides the function definitions for the ShallowWaterProblemParameters
 *        class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] problem_name_     name of problem
 * \param[in] is_steady_state_  flag that problem is to be run in steady-state
 */
template <int dim>
ShallowWaterProblemParameters<dim>::ShallowWaterProblemParameters(
  const std::string & problem_name_, const bool & is_steady_state_)
  : ProblemParameters<dim>(problem_name_, dim + 1, is_steady_state_)
{
}

/**
 * \brief Declares parameters
 */
template <int dim>
void ShallowWaterProblemParameters<dim>::declare_derived_parameters()
{
  // bathymetry
  this->parameter_handler.enter_subsection("bathymetry");
  {
    this->parameter_handler.declare_entry(
      "bathymetry function", "0", Patterns::Anything(), "Bathymetry function");
  }
  this->parameter_handler.leave_subsection();

  // boundary conditions
  this->parameter_handler.enter_subsection("boundary conditions");
  {
    this->parameter_handler.declare_entry("dirichlet function height",
                                          "0",
                                          Patterns::Anything(),
                                          "Dirichlet function for height");
    this->parameter_handler.declare_entry("dirichlet function momentumx",
                                          "0",
                                          Patterns::Anything(),
                                          "Dirichlet function for x-momentum");
    this->parameter_handler.declare_entry("dirichlet function momentumy",
                                          "0",
                                          Patterns::Anything(),
                                          "Dirichlet function for y-momentum");
  }
  this->parameter_handler.leave_subsection();

  // initial conditions
  this->parameter_handler.enter_subsection("initial conditions");
  {
    this->parameter_handler.declare_entry("initial conditions height",
                                          "0",
                                          Patterns::Anything(),
                                          "Initial conditions for height");
    this->parameter_handler.declare_entry("initial conditions momentumx",
                                          "0",
                                          Patterns::Anything(),
                                          "Initial conditions for x-momentum");
    this->parameter_handler.declare_entry("initial conditions momentumy",
                                          "0",
                                          Patterns::Anything(),
                                          "Initial conditions for y-momentum");
  }
  this->parameter_handler.leave_subsection();

  // exact solution
  this->parameter_handler.enter_subsection("exact solution");
  {
    this->parameter_handler.declare_entry("exact solution height",
                                          "1",
                                          Patterns::Anything(),
                                          "Exact solution for height");
    this->parameter_handler.declare_entry("exact solution momentumx",
                                          "1",
                                          Patterns::Anything(),
                                          "Exact solution for x-momentum");
    this->parameter_handler.declare_entry("exact solution momentumy",
                                          "1",
                                          Patterns::Anything(),
                                          "Exact solution for y-momentum");
  }
  this->parameter_handler.leave_subsection();

  // constants
  this->parameter_handler.enter_subsection("constants");
  {
    this->parameter_handler.declare_entry(
      "gravity", "1.0", Patterns::Double(), "Acceleration due to gravity");
    this->parameter_handler.declare_entry(
      "x_interface", "0.0", Patterns::Double(), "x-position of interface");
    this->parameter_handler.declare_entry(
      "h_left", "1.0", Patterns::Double(), "Left height value");
    this->parameter_handler.declare_entry(
      "h_right", "1.0", Patterns::Double(), "Right height value");
    this->parameter_handler.declare_entry(
      "h_unperturbed", "1.0", Patterns::Double(), "Unperturbed height value");
    this->parameter_handler.declare_entry(
      "h_perturbed", "1.0", Patterns::Double(), "Perturbed height value");
    this->parameter_handler.declare_entry(
      "u_left", "0.0", Patterns::Double(), "Left x-velocity value");
    this->parameter_handler.declare_entry(
      "u_right", "0.0", Patterns::Double(), "Right x-velocity value");
    this->parameter_handler.declare_entry("bump_x_center",
                                          "0.0",
                                          Patterns::Double(),
                                          "Center of bump in x-direction");
    this->parameter_handler.declare_entry("bump_y_center",
                                          "0.0",
                                          Patterns::Double(),
                                          "Center of bump in y-direction");
    this->parameter_handler.declare_entry(
      "bump_x_width", "0.0", Patterns::Double(), "Width of bump in x-direction");
    this->parameter_handler.declare_entry(
      "bump_y_width", "0.0", Patterns::Double(), "Width of bump in y-direction");
    this->parameter_handler.declare_entry(
      "bump_height", "0.0", Patterns::Double(), "Height of bump");
    this->parameter_handler.declare_entry(
      "perturbation_x_center",
      "0.0",
      Patterns::Double(),
      "Center of perturbation in x-direction");
    this->parameter_handler.declare_entry(
      "perturbation_y_center",
      "0.0",
      Patterns::Double(),
      "Center of perturbation in y-direction");
    this->parameter_handler.declare_entry("perturbation_x_width",
                                          "0.0",
                                          Patterns::Double(),
                                          "Width of perturbation in x-direction");
    this->parameter_handler.declare_entry("perturbation_y_width",
                                          "0.0",
                                          Patterns::Double(),
                                          "Width of perturbation in y-direction");
  }
  this->parameter_handler.leave_subsection();
}

/**
 * \brief Gets parameters from parameter handler.
 */
template <int dim>
void ShallowWaterProblemParameters<dim>::get_derived_parameters()
{
  // bathymetry
  this->parameter_handler.enter_subsection("bathymetry");
  {
    bathymetry_function_string =
      this->parameter_handler.get("bathymetry function");
  }
  this->parameter_handler.leave_subsection();

  // boundary conditions
  this->parameter_handler.enter_subsection("boundary conditions");
  {
    dirichlet_function_height =
      this->parameter_handler.get("dirichlet function height");
    dirichlet_function_momentumx =
      this->parameter_handler.get("dirichlet function momentumx");
    dirichlet_function_momentumy =
      this->parameter_handler.get("dirichlet function momentumy");
  }
  this->parameter_handler.leave_subsection();

  // initial conditions
  this->parameter_handler.enter_subsection("initial conditions");
  {
    initial_conditions_height =
      this->parameter_handler.get("initial conditions height");
    initial_conditions_momentumx =
      this->parameter_handler.get("initial conditions momentumx");
    initial_conditions_momentumy =
      this->parameter_handler.get("initial conditions momentumy");
  }
  this->parameter_handler.leave_subsection();

  // exact solution
  this->parameter_handler.enter_subsection("exact solution");
  {
    exact_solution_height = this->parameter_handler.get("exact solution height");
    exact_solution_momentumx =
      this->parameter_handler.get("exact solution momentumx");
    exact_solution_momentumy =
      this->parameter_handler.get("exact solution momentumy");
  }
  this->parameter_handler.leave_subsection();

  // constants
  this->parameter_handler.enter_subsection("constants");
  {
    gravity = this->parameter_handler.get_double("gravity");
    x_interface = this->parameter_handler.get_double("x_interface");
    h_left = this->parameter_handler.get_double("h_left");
    h_right = this->parameter_handler.get_double("h_right");
    h_unperturbed = this->parameter_handler.get_double("h_unperturbed");
    h_perturbed = this->parameter_handler.get_double("h_perturbed");
    u_left = this->parameter_handler.get_double("u_left");
    u_right = this->parameter_handler.get_double("u_right");
    bump_x_center = this->parameter_handler.get_double("bump_x_center");
    bump_y_center = this->parameter_handler.get_double("bump_y_center");
    bump_x_width = this->parameter_handler.get_double("bump_x_width");
    bump_y_width = this->parameter_handler.get_double("bump_y_width");
    bump_height = this->parameter_handler.get_double("bump_height");
    perturbation_x_center =
      this->parameter_handler.get_double("perturbation_x_center");
    perturbation_y_center =
      this->parameter_handler.get_double("perturbation_y_center");
    perturbation_x_width =
      this->parameter_handler.get_double("perturbation_x_width");
    perturbation_y_width =
      this->parameter_handler.get_double("perturbation_y_width");
  }
  this->parameter_handler.leave_subsection();
}

/**
 * \brief Processes derived parameters.
 *
 * \param[in] triangulation    triangulation
 * \param[in] fe               finite element system
 * \param[in] face_quadrature  face quadrature
 */
template <int dim>
void ShallowWaterProblemParameters<dim>::process_derived_parameters(
  Triangulation<dim> &,
  const FESystem<dim> & fe,
  const QGauss<dim - 1> & face_quadrature)
{
  // constants for function parsers
  this->constants["gravity"] = gravity;
  this->constants["x_interface"] = x_interface;
  this->constants["h_left"] = h_left;
  this->constants["h_right"] = h_right;
  this->constants["h_unperturbed"] = h_unperturbed;
  this->constants["h_perturbed"] = h_perturbed;
  this->constants["u_left"] = u_left;
  this->constants["u_right"] = u_right;
  this->constants["bump_height"] = bump_height;
  this->constants["bump_x_center"] = bump_x_center;
  this->constants["bump_y_center"] = bump_y_center;
  this->constants["bump_x_width"] = bump_x_width;
  this->constants["bump_y_width"] = bump_y_width;
  this->constants["bump_left"] = bump_x_center - 0.5 * bump_x_width;
  this->constants["bump_right"] = bump_x_center + 0.5 * bump_x_width;
  this->constants["perturbation_x_center"] = perturbation_x_center;
  this->constants["perturbation_y_center"] = perturbation_y_center;
  this->constants["perturbation_x_width"] = perturbation_x_width;
  this->constants["perturbation_y_width"] = perturbation_y_width;

  // store Dirichlet BC strings
  this->dirichlet_function_strings.resize(this->n_components);
  this->dirichlet_function_strings[0] = dirichlet_function_height;
  this->dirichlet_function_strings[1] = dirichlet_function_momentumx;
  if (dim == 2)
    this->dirichlet_function_strings[2] = dirichlet_function_momentumy;

  // store initial condition strings
  this->initial_conditions_strings.resize(this->n_components);
  this->initial_conditions_strings[0] = initial_conditions_height;
  this->initial_conditions_strings[1] = initial_conditions_momentumx;
  if (dim == 2)
    this->initial_conditions_strings[2] = initial_conditions_momentumy;

  // exact solution
  if (this->has_exact_solution)
  {
    if (this->exact_solution_type == "riemann")
    {
      // create and initialize Riemann solver for exact solution
      std::shared_ptr<ShallowWaterRiemannSolver<dim>>
        exact_solution_function_derived =
          std::make_shared<ShallowWaterRiemannSolver<dim>>(
            h_left, u_left, h_right, u_right, gravity, x_interface);
      this->exact_solution_function = exact_solution_function_derived;
    }
    else // assumed to be function parser
    {
      // store exact solution strings
      this->exact_solution_strings.resize(this->n_components);
      this->exact_solution_strings[0] = exact_solution_height;
      this->exact_solution_strings[1] = exact_solution_momentumx;
      if (dim == 2)
        this->exact_solution_strings[2] = exact_solution_momentumy;
    }
  }

  // boundary conditions
  if (this->boundary_conditions_type == "none")
  {
    std::shared_ptr<ShallowWaterNoBC<dim>> derived_boundary_conditions =
      std::make_shared<ShallowWaterNoBC<dim>>(fe, face_quadrature, gravity);
    this->boundary_conditions = derived_boundary_conditions;
  }
  else if (this->boundary_conditions_type == "2d_dam_break")
  {
    std::shared_ptr<ShallowWater2DDamBreakBC<dim>> derived_boundary_conditions =
      std::make_shared<ShallowWater2DDamBreakBC<dim>>(
        fe, face_quadrature, gravity);
    this->boundary_conditions = derived_boundary_conditions;
  }
  else if (this->boundary_conditions_type == "characteristic_open")
  {
    std::shared_ptr<ShallowWaterSubcriticalOpenBC1D<dim>>
      derived_boundary_conditions =
        std::make_shared<ShallowWaterSubcriticalOpenBC1D<dim>>(
          fe, face_quadrature, gravity, h_unperturbed, h_unperturbed);
    this->boundary_conditions = derived_boundary_conditions;
  }
  else if (this->boundary_conditions_type == "wall")
  {
    std::shared_ptr<ShallowWaterWallBC<dim>> derived_boundary_conditions =
      std::make_shared<ShallowWaterWallBC<dim>>(fe, face_quadrature, gravity);
    this->boundary_conditions = derived_boundary_conditions;
  }
  else if (this->boundary_conditions_type == "characteristic_wall")
  {
    std::shared_ptr<ShallowWaterSubcriticalWallBC1D<dim>>
      derived_boundary_conditions =
        std::make_shared<ShallowWaterSubcriticalWallBC1D<dim>>(
          fe, face_quadrature, gravity);
    this->boundary_conditions = derived_boundary_conditions;
  }

  // bathymetry function
  std::shared_ptr<FunctionParser<dim>> bathymetry_function_derived =
    std::make_shared<FunctionParser<dim>>(1);
  std::map<std::string, double> constants;
  bathymetry_function_derived->initialize(
    FunctionParser<dim>::default_variable_names(),
    bathymetry_function_string,
    this->constants,
    false);
  bathymetry_function = bathymetry_function_derived;
}
