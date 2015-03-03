/** \brief constructor
 */
TransportParameters::TransportParameters() :
      problem_id(1),
      degree(1),
      time_integrator_option(1),
      refinement_option(1),
      time_refinement_factor(0.5),
      use_adaptive_mesh_refinement(false),
      n_refinement_cycles(5),
      initial_refinement_level(2),
      linear_solver_option(1),
      scheme_option(0),
      low_order_diffusion_option(1),
      entropy_string("0.5*u*u"),
      entropy_derivative_string("u"),
      entropy_residual_coefficient(1.0),
      jump_coefficient(1.0),
      output_meshes(false),
      is_steady_state(true),
      end_time(1.0),
      time_step_size(0.001),
      output_exact_solution(false),
      output_initial_solution(false),
      output_DMP_bounds(false),
      CFL_limit(0.5),
      do_not_limit(false),
      n_quadrature_points(3),
      save_convergence_results(false)
{
}

/** \brief defines all of the input parameters
 */
void TransportParameters::declare_parameters(ParameterHandler &prm)
{
   prm.declare_entry("Problem ID", "1", Patterns::Integer(),
         "Problem ID");
   prm.declare_entry("Finite element degree", "1", Patterns::Integer(),
         "Polynomial degree of finite elements");
   prm.declare_entry("Time integrator option", "1", Patterns::Integer(),
         "Choice of time integrator");
   prm.declare_entry("Refinement option", "1", Patterns::Integer(),
         "Refinement option (space and/or time)");
   prm.declare_entry("Time refinement factor", "0.5", Patterns::Double(),
         "Factor by which to reduce dt if refining time step size");
   prm.declare_entry("Use adaptive mesh refinement", "false", Patterns::Bool(),
         "Option to use adaptive mesh refinement instead of uniform");
   prm.declare_entry("Number of refinement cycles", "5", Patterns::Integer(),
         "Number of mesh refinement cycles");
   prm.declare_entry("Initial refinement level", "2", Patterns::Integer(),
         "Number of refinements for first mesh refinement cycle");
   prm.declare_entry("Linear solver option", "1", Patterns::Integer(),
         "Option for linear solver");
   prm.declare_entry("Scheme option", "0", Patterns::Integer(),
         "Option for scheme to be used");
   prm.declare_entry("Low order diffusion option", "1", Patterns::Integer(),
         "Option for low-order diffusion to be used");
   prm.declare_entry("Entropy string", "0.5*u*u", Patterns::Anything(),
         "String for entropy function");
   prm.declare_entry("Entropy derivative string", "0.5*u*u", Patterns::Anything(),
         "String for entropy derivative function");
   prm.declare_entry("Entropy residual coefficient", "1.0", Patterns::Double(),
         "Coefficient for the entropy viscosity");
   prm.declare_entry("Jump coefficient", "1.0", Patterns::Double(),
         "Coefficient for jumps used with entropy viscosity");
   prm.declare_entry("Output mesh", "false", Patterns::Bool(),
         "Option to output meshes as .eps files");
   prm.declare_entry("Is steady state", "true", Patterns::Bool(),
         "Boolean flag for steady-state vs. transient");
   prm.declare_entry("End time", "1.0", Patterns::Double(),
         "End time if transient problem is run");
   prm.declare_entry("Time step size", "0.01", Patterns::Double(),
         "Time step size if transient problem is run");
   prm.declare_entry("Output exact solution", "false", Patterns::Bool(),
         "Option to output exact solution");
   prm.declare_entry("Output initial solution", "false", Patterns::Bool(),
         "Option to output initial solution");
   prm.declare_entry("Output DMP bounds", "false", Patterns::Bool(),
         "Option to output DMP bounds");
   prm.declare_entry("CFL limit", "0.5", Patterns::Double(),
         "Upper bound for the CFL number");
   prm.declare_entry("Do not limit", "false", Patterns::Bool(),
         "Option to choose not to limit when constructing high-order solution");
   prm.declare_entry("Number of quadrature points", "3", Patterns::Integer(),
         "Number of quadrature points to use in formula");
   prm.declare_entry("Save convergence results", "false", Patterns::Bool(),
         "Option save convegence results");
}

/** \brief get the input parameters
 */
void TransportParameters::get_parameters(ParameterHandler &prm) {
   problem_id = prm.get_integer("Problem ID");
   degree = prm.get_integer("Finite element degree");
   time_integrator_option = prm.get_integer("Time integrator option");
   refinement_option = prm.get_integer("Refinement option");
   time_refinement_factor = prm.get_double("Time refinement factor");
   use_adaptive_mesh_refinement = prm.get_bool("Use adaptive mesh refinement");
   n_refinement_cycles = prm.get_integer("Number of refinement cycles");
   initial_refinement_level = prm.get_integer("Initial refinement level");
   linear_solver_option = prm.get_integer("Linear solver option");
   scheme_option = prm.get_integer("Scheme option");
   low_order_diffusion_option = prm.get_integer("Low order diffusion option");
   entropy_string = prm.get("Entropy string");
   entropy_derivative_string = prm.get("Entropy derivative string");
   entropy_residual_coefficient = prm.get_double("Entropy residual coefficient");
   jump_coefficient = prm.get_double("Jump coefficient");
   output_meshes = prm.get_bool("Output mesh");
   is_steady_state = prm.get_bool("Is steady state");
   end_time = prm.get_double("End time");
   time_step_size = prm.get_double("Time step size");
   output_exact_solution = prm.get_bool("Output exact solution");
   output_initial_solution = prm.get_bool("Output initial solution");
   output_DMP_bounds = prm.get_bool("Output DMP bounds");
   CFL_limit = prm.get_double("CFL limit");
   do_not_limit = prm.get_bool("Do not limit");
   n_quadrature_points = prm.get_integer("Number of quadrature points");
   save_convergence_results = prm.get_bool("Save convergence results");
}
