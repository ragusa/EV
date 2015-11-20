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
 * \brief Declares common problem parameters.
 *
 * \param[out] parameter_handler parameter handler
 */
template <int dim>
void ProblemParameters<dim>::declare_common_parameters(
  ParameterHandler & parameter_handler)
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
      "x width", "1.0", Patterns::Double(), "Width of domain in x-direction");
    parameter_handler.declare_entry(
      "y width", "1.0", Patterns::Double(), "Width of domain in y-direction");
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
 * \brief Gets common problem parameters from parameter handler.
 *
 * \param[in] parameter_handler parameter handler
 */
template <int dim>
void ProblemParameters<dim>::get_common_parameters(
  ParameterHandler & parameter_handler)
{
  // dimension
  parameter_handler.enter_subsection("dimension");
  {
    valid_in_1d = parameter_handler.get_bool("valid in 1d");
    valid_in_2d = parameter_handler.get_bool("valid in 2d");
  }
  parameter_handler.leave_subsection();

  // domain
  parameter_handler.enter_subsection("domain");
  {
    domain_shape = parameter_handler.get("domain shape");
    x_start = parameter_handler.get_double("x start");
    y_start = parameter_handler.get_double("y start");
    x_width = parameter_handler.get_double("x width");
    y_width = parameter_handler.get_double("y width");
  }
  parameter_handler.leave_subsection();

  // boundary conditions
  parameter_handler.enter_subsection("boundary conditions");
  {
    boundary_conditions_type = parameter_handler.get("boundary conditions type");
    use_exact_solution_as_dirichlet_bc =
      parameter_handler.get_bool("use exact solution as dirichlet bc");
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
