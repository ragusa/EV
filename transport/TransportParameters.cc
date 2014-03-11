/** \brief constructor
 */
//template<int dim>
TransportParameters::TransportParameters() :
      degree(1),
      n_energy_groups(1),
      n_directions(1),
      use_adaptive_mesh_refinement(false),
      n_refinement_cycles(5),
      initial_refinement_level(2),
      solver_option(1),
      preconditioner_option(1),
      source_option(1),
      source_value(0.0e0),
      total_cross_section_option(1),
      total_cross_section_value(1.0e0),
      incoming_flux(1.0e0),
      viscosity_type(0),
      max_viscosity_coefficient(5.0e-1),
      entropy_viscosity_coefficient(1.0e-1),
      max_nonlinear_iterations(10),
      relative_difference_tolerance(1.0e-6),
      exact_solution_id(0),
      output_meshes(false)
{
}

/** \brief defines all of the input parameters
 */
//template<int dim>
void TransportParameters::declare_parameters(ParameterHandler &prm)
{
   prm.declare_entry("Finite element degree", "1", Patterns::Integer(),
         "Polynomial degree of finite elements");
   prm.declare_entry("Number of energy groups", "1", Patterns::Integer(),
         "Number of neutron energy groups");
   prm.declare_entry("Number of directions", "1", Patterns::Integer(),
         "Number of transport directions");
   prm.declare_entry("Use adaptive mesh refinement", "false",
         Patterns::Bool(),
         "Option to use adaptive mesh refinement instead of uniform");
   prm.declare_entry("Number of refinement cycles", "5", Patterns::Integer(),
         "Number of mesh refinement cycles");
   prm.declare_entry("Initial refinement level", "2", Patterns::Integer(),
         "Number of refinements for first mesh refinement cycle");
   prm.declare_entry("Solver option", "1", Patterns::Integer(),
         "Option for linear solver");
   prm.declare_entry("Preconditioner option", "1", Patterns::Integer(),
         "Option for preconditioner for linear solver");
   prm.declare_entry("Source option", "1", Patterns::Integer(),
         "Option for source definition");
   prm.declare_entry("Source value", "0.0e0", Patterns::Double(),
         "Value of extraneous source term");
   prm.declare_entry("Total cross section option", "1", Patterns::Integer(),
         "Option for total cross section definition");
   prm.declare_entry("Total cross section value", "1.0e0", Patterns::Double(),
         "Value of total cross section");
   prm.declare_entry("Incoming flux", "1.0e0", Patterns::Double(),
         "Value for flux on incoming boundary");
   prm.declare_entry("Viscosity type", "0", Patterns::Integer(),
         "Option for viscosity type: none, first-order, or entropy");
   prm.declare_entry("Max viscosity coefficient", "5.0e-1", Patterns::Double(),
         "Coefficient for the first-order viscosity");
   prm.declare_entry("Entropy viscosity coefficient", "1.0e-1", Patterns::Double(),
         "Coefficient for the entropy viscosity");
   prm.declare_entry("Maximum number of nonlinear iterations", "10", Patterns::Integer(),
         "Maximum number of nonlinear iterations to use when using entropy viscosity");
   prm.declare_entry("Relative difference tolerance", "1.0e-6", Patterns::Double(),
         "Relative difference tolerance for nonlinear convergence");
   prm.declare_entry("Exact solution ID", "0", Patterns::Integer(),
         "ID of exact solution to use when evaluating error");
   prm.declare_entry("Output mesh", "false", Patterns::Bool(),
         "Option to output meshes as .eps files");
}

/** \brief get the input parameters
 */
//template<int dim>
void TransportParameters::get_parameters(ParameterHandler &prm) {
   degree = prm.get_integer("Finite element degree");
   n_energy_groups = prm.get_integer("Number of energy groups");
   n_directions = prm.get_integer("Number of directions");
   use_adaptive_mesh_refinement = prm.get_bool("Use adaptive mesh refinement");
   n_refinement_cycles = prm.get_integer("Number of refinement cycles");
   initial_refinement_level = prm.get_integer("Initial refinement level");
   solver_option = prm.get_integer("Solver option");
   preconditioner_option = prm.get_integer("Preconditioner option");
   source_option = prm.get_integer("Source option");
   source_value = prm.get_double("Source value");
   total_cross_section_option = prm.get_integer("Total cross section option");
   total_cross_section_value = prm.get_double("Total cross section value");
   incoming_flux = prm.get_double("Incoming flux");
   viscosity_type = prm.get_integer("Viscosity type");
   max_viscosity_coefficient = prm.get_double(
         "Max viscosity coefficient");
   entropy_viscosity_coefficient = prm.get_double(
         "Entropy viscosity coefficient");
   max_nonlinear_iterations = prm.get_integer("Maximum number of nonlinear iterations");
   relative_difference_tolerance = prm.get_double("Relative difference tolerance");
   exact_solution_id = prm.get_integer("Exact solution ID");
   output_meshes = prm.get_bool("Output mesh");
}
