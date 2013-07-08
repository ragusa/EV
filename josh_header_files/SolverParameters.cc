#include <deal.II/base/parameter_handler.h>

#include "SolverParameters.h"
	
// ----------------------------------
// --- DECLARE SOLVER PARAMETERS
// ----------------------------------
void SolverParameters::declare_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("nonlinear solver");
      {
        prm.declare_entry("nonlinear verbosity", "quiet",
                          Patterns::Selection("quiet|verbose"),
                          "State whether output from nonlinear solver runs should be printed. "
                          "Choices are <quiet|verbose>.");
        prm.declare_entry("nonlinear method", "newton",
                          Patterns::Selection("newton"),
                          "The kind of nonlinear solver for the linear system. "
                          "Choices are <newton>.");
        prm.declare_entry("nonlinear absolute tolerance", "1e-10",
                          Patterns::Double(),
                          "Nonlinear absolute tolerance");
        prm.declare_entry("nonlinear relative tolerance", "1e-10",
                          Patterns::Double(),
                          "Nonlinear relative tolerance");
        prm.declare_entry("max nonlinear iterations", "300",
                          Patterns::Integer(),
                          "Maximum nonlinear iterations");
        prm.declare_entry("damping", "1.0",
                          Patterns::Double(),
                          "damping");
      }
    prm.leave_subsection();
	
    prm.enter_subsection("linear solver");
      {
        prm.declare_entry("linear verbosity", "quiet",
                          Patterns::Selection("quiet|verbose"),
                          "State whether output from linear solver runs should be printed. "
                          "Choices are <quiet|verbose>.");
        prm.declare_entry("linear method", "gmres",
                          Patterns::Selection("gmres|direct"),
                          "The kind of linear solver for the linear system. "
                          "Choices are <gmres|direct>.");
        prm.declare_entry("linear absolute tolerance", "1e-10",
                          Patterns::Double(),
                          "Linear absolute tolerance");
        prm.declare_entry("linear relative tolerance", "1e-10",
                          Patterns::Double(),
                          "Linear relative tolerance");
        prm.declare_entry("max linear iterations", "300",
                          Patterns::Integer(),
                          "Maximum linear solver iterations");
/*
        prm.declare_entry("ilut fill", "2",
                          Patterns::Double(),
                          "Ilut preconditioner fill");
        prm.declare_entry("ilut absolute tolerance", "1e-9",
                          Patterns::Double(),
                          "Ilut preconditioner tolerance");
        prm.declare_entry("ilut relative tolerance", "1.1",
                          Patterns::Double(),
                          "Ilut relative tolerance");
        prm.declare_entry("ilut drop tolerance", "1e-10",
                          Patterns::Double(),
                          "Ilut drop tolerance");
*/
      }
    prm.leave_subsection();
}

// ----------------------------------
// --- PARSE SOLVER PARAMETERS
// ----------------------------------
void SolverParameters::parse_parameters (ParameterHandler &prm)
{
    prm.enter_subsection("nonlinear solver");
      {
        const std::string verbosity = prm.get("nonlinear verbosity");
        if (verbosity == "verbose") nonlinear_verbosity = verbose;
        if (verbosity == "quiet")   nonlinear_verbosity = quiet;
          
        const std::string solver = prm.get("nonlinear method");
        if (solver == "newton")     nonlinear_solver = newton;
          
        nonlinear_atol         = prm.get_double("nonlinear absolute tolerance");
        nonlinear_rtol         = prm.get_double("nonlinear relative tolerance");
        max_nonlin_iterations  = prm.get_integer("max nonlinear iterations");
        damping                = prm.get_double("damping");
      }
    prm.leave_subsection();
	
    prm.enter_subsection("linear solver");
      {
        const std::string verbosity = prm.get("linear verbosity");
        if (verbosity == "verbose") linear_verbosity = verbose;
        if (verbosity == "quiet")   linear_verbosity = quiet;
          
        const std::string solver = prm.get("linear method");
        if (solver == "direct")     linear_solver = direct;
        else if (solver == "gmres") linear_solver = gmres;
          
        linear_atol     = prm.get_double("linear absolute tolerance");
        linear_rtol     = prm.get_double("linear relative tolerance");
        max_linear_iterations  = prm.get_integer("max linear iterations");
/*
        ilut_fill       = prm.get_double("ilut fill");
        ilut_atol       = prm.get_double("ilut absolute tolerance");
        ilut_rtol       = prm.get_double("ilut relative tolerance");
        ilut_drop       = prm.get_double("ilut drop tolerance");
*/
      }
    prm.leave_subsection();
}