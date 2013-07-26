/** \fn declare_parameters
 *  \brief Declares the parameters that will be in input file.
 *
 *  This function declares all of the input parameters required
 *  for the ConservationLaw class.
 */
template <int dim>
void ConservationLawParameters<dim>::declare_parameters (ParameterHandler &prm)
{
	// time parameters
	prm.enter_subsection("time");
	{
		prm.declare_entry("final time", "1.0",
				Patterns::Double(),
				"final time value");
		prm.declare_entry("time step size", "1e-3",
				Patterns::Double(),
				"time step size");
	}
        prm.leave_subsection();

    // temporal integrator
    prm.enter_subsection("temporal integrator");
    {
        prm.declare_entry("temporal integrator", "erk",
                          Patterns::Selection("erk"),
                          "The method used for advancing time. "
                          "Choices are <erk>.");
        prm.declare_entry("erk stages","1",
                          Patterns::Integer(),
                          "Number of stages for explicit Runge-Kutta.");
    }
    prm.leave_subsection();

    // nonlinear solver parameters
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
	
    // linear solver parameters
    prm.enter_subsection("linear solver");
      {
        prm.declare_entry("linear verbosity", "quiet",
                          Patterns::Selection("quiet|verbose"),
                          "State whether output from linear solver runs should be printed. "
                          "Choices are <quiet|verbose>.");
        prm.declare_entry("linear method", "gmres",
                          Patterns::Selection("gmres|direct|bicgstab"),
                          "The kind of linear solver for the linear system. "
                          "Choices are <gmres|direct|bicgstab>.");
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

/*
   // initial conditions
   prm.enter_subsection("initial conditions");
   {
      for (int c = 0; c < n_components; ++c)
         prm.declare_entry("initial conditions " + Utilities::int_to_string(c),
         "0.0",
         Patterns::Anything(),
         "initial conditions for component computed from x,y,z");
   }
   prm.leave_subsection();
*/
    prm.enter_subsection("output");
    {
    	prm.declare_entry("output period", "1",
    			Patterns::Integer(),
    			"Period of time steps for outputting the solution, e.g.,"
    			" 1 would output every time step,"
    			" and 2 would output every other time step, etc.");
    }
    prm.leave_subsection();
}

/** \fn get_parameters
 *  \brief Gets the parameters from the parameter handler.
 *
 *  This function takes the input parameters from the parameter
 *  handler into the member variables.
 */
template <int dim>
void ConservationLawParameters<dim>::get_parameters (ParameterHandler &prm)
{
	// time parameters
	prm.enter_subsection("time");
	{
		final_time = prm.get_double("final time");
		time_step_size = prm.get_double("time step size");
	}
	prm.leave_subsection();

    // temporal integrator
    prm.enter_subsection("temporal integrator");
    {
        const std::string temporal_choice = prm.get("temporal integrator");
        if (temporal_choice == "erk")
           temporal_integrator = erk;
        else
           Assert(false,ExcNotImplemented());

        erk_nstages = prm.get_integer("erk stages");
    }
    prm.leave_subsection();

    // nonlinear solver parameters
    prm.enter_subsection("nonlinear solver");
      {
        const std::string verbosity = prm.get("nonlinear verbosity");
        if (verbosity == "verbose") nonlinear_verbosity = verbose;
        if (verbosity == "quiet")   nonlinear_verbosity = quiet;
          
        const std::string solver = prm.get("nonlinear method");
        if (solver == "newton")     nonlinear_solver = newton;
          
        nonlinear_atol            = prm.get_double("nonlinear absolute tolerance");
        nonlinear_rtol            = prm.get_double("nonlinear relative tolerance");
        max_nonlinear_iterations  = prm.get_integer("max nonlinear iterations");
        damping                   = prm.get_double("damping");
      }
    prm.leave_subsection();
	
    // linear solver parameters
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

/*
   // initial conditions
   prm.enter_subsection("initial conditions");
   {
      std::vector<std::string> expressions (n_components,"0.0");
      for (int c = 0; c < n_components; c++)
          expressions[c] = prm.get("initial conditions " + Utilities::int_to_string(c));
      initial_conditions.initialize (FunctionParser<dim>::default_variable_names(),
                                     expressions,
                                     std::map<std::string, double>());
   }
   prm.leave_subsection();
*/
}
