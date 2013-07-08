#include <deal.II/base/parameter_handler.h>

class SolverParameters
{
  public:
    SolverParameters()
	~SolverParameters();
  
    enum NonlinearSolverType { newton };
    NonlinearSolverType nonlinear_solver;
	
    enum LinearSolverType { gmres, direct };
    LinearSolverType linear_solver;
    
    enum Verbosity { quiet, verbose };
    Verbosity nonlinear_verbosity;
    Verbosity linear_verbosity;
      
    double linear_atol;
    double linear_rtol;
    int max_linear_iterations;
/*
    double ilut_fill;
    double ilut_atol;
    double ilut_rtol;
    double ilut_drop;
*/
    
    double nonlinear_atol;
    double nonlinear_rtol;
    int max_nonlin_iterations;
    double damping;
	
	void declare_parameters (ParameterHandler &prm);
    void parse_parameters (ParameterHandler &prm);
};