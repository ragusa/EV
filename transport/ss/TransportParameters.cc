/** \brief constructor
 */
template <int dim>
TransportParameters<dim>::TransportParameters()
  : problem_id(1),
    degree(1),
    use_adaptive_refinement(false),
    n_refinement_cycles(5),
    initial_refinement_level(2),
    linear_solver_option(1),
    output_solution(false),
    n_quadrature_points(3)
{
}

/** \brief defines all of the input parameters
 */
template <int dim>
void TransportParameters<dim>::declare_parameters(ParameterHandler & prm)
{
  prm.declare_entry("Problem ID", "1", Patterns::Integer(), "Problem ID");
  prm.declare_entry("Finite element degree",
                    "1",
                    Patterns::Integer(),
                    "Polynomial degree of finite elements");
  prm.declare_entry("Use adaptive refinement",
                    "false",
                    Patterns::Bool(),
                    "Option to use adaptive mesh refinement instead of uniform");
  prm.declare_entry("Number of refinement cycles",
                    "5",
                    Patterns::Integer(),
                    "Number of mesh refinement cycles");
  prm.declare_entry("Initial refinement level",
                    "2",
                    Patterns::Integer(),
                    "Number of refinements for first mesh refinement cycle");
  prm.declare_entry(
    "Linear solver option", "1", Patterns::Integer(), "Option for linear solver");
  prm.declare_entry(
    "Output solution", "false", Patterns::Bool(), "Option to output solution");
  prm.declare_entry("Number of quadrature points",
                    "3",
                    Patterns::Integer(),
                    "Number of quadrature points to use in formula");
}

/** \brief get the input parameters
 */
template <int dim>
void TransportParameters<dim>::get_parameters(ParameterHandler & prm)
{
  problem_id = prm.get_integer("Problem ID");
  degree = prm.get_integer("Finite element degree");
  use_adaptive_refinement = prm.get_bool("Use adaptive refinement");
  n_refinement_cycles = prm.get_integer("Number of refinement cycles");
  initial_refinement_level = prm.get_integer("Initial refinement level");
  linear_solver_option = prm.get_integer("Linear solver option");
  output_solution = prm.get_bool("Output solution");
  n_quadrature_points = prm.get_integer("Number of quadrature points");
}
