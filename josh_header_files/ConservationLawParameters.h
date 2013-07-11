#ifndef ConservationLawParameters_h
#define ConservationLawParameters_h

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>

using namespace dealii;

template<int dim>
class ConservationLawParameters
{
  public:
    ConservationLawParameters(const int &n_components);
//    ~ConservationLawParameters();
  
    enum NonlinearSolverType { newton };
    NonlinearSolverType nonlinear_solver;
	
    enum LinearSolverType { gmres, direct };
    LinearSolverType linear_solver;
    
    enum Verbosity { quiet, verbose };
    Verbosity nonlinear_verbosity;
    Verbosity linear_verbosity;
      
    int n_components;
    FunctionParser<dim> initial_conditions;

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
    void get_parameters (ParameterHandler &prm);
};

#include "ConservationLawParameters.cc"

#endif