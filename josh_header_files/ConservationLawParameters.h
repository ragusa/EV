/** \file ConservationLawParameters.h
 *  \brief Provides the header for the ConservationLawParameters class.
 */
#ifndef ConservationLawParameters_h
#define ConservationLawParameters_h

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/** \class ConservationLawParameters
 *  \brief Declares and retrieves input parameters related to solving
 *  for a conservation law.
 */
template<int dim>
class ConservationLawParameters
{
  public:

    int n_components;

    enum TimeStepSizeMethod { constant, cfl_condition };
    TimeStepSizeMethod time_step_size_method;
    double final_time;
    double time_step_size;
    double cfl;

    unsigned int output_period;
  
    enum TemporalIntegrator { erk };
    TemporalIntegrator temporal_integrator;
    int erk_nstages; 

    enum NonlinearSolverType { newton };
    NonlinearSolverType nonlinear_solver;

    enum LinearSolverType { direct, gmres };
    LinearSolverType linear_solver;
    LinearSolverType mass_matrix_linear_solver;
    
    enum Verbosity { quiet, verbose };
    Verbosity nonlinear_verbosity;
    Verbosity linear_verbosity;
      
    double linear_atol;
    double linear_rtol;
    int max_linear_iterations;
    
    double nonlinear_atol;
    double nonlinear_rtol;
    int max_nonlinear_iterations;
    double damping;

    static void declare_conservation_law_parameters (ParameterHandler &prm);
    void get_conservation_law_parameters (ParameterHandler &prm);

    int degree;
};

#include "ConservationLawParameters.cc"

#endif
