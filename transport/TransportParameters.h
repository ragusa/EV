#ifndef TransportParameters_cc
#define TransportParameters_cc

#include <deal.II/base/parameter_handler.h>
#include "SSPRKTimeIntegrator.h"
#include "RefinementHandler.h"
#include "EntropyViscosity.h"

using namespace dealii;

/** \brief parameters class for defining input parameters
*/
template <int dim>
class TransportParameters {
 public:
    TransportParameters();
    static void declare_parameters(ParameterHandler &prm);
    void get_parameters(ParameterHandler &prm);

    unsigned int problem_id; // problem ID

    unsigned int degree; // polynomial degree of finite elements
    typename SSPRKTimeIntegrator<dim>::SSPRKMethod time_integrator_option; // SSPRK method
    typename RefinementHandler<dim>::RefinementMode refinement_mode; // refinement mode (space or time)
    double time_refinement_factor; // reduction factor for time refinement
    bool use_adaptive_refinement; // option to use adaptive mesh refinement
    unsigned int n_refinement_cycles; // number of refinement cycles
    unsigned int initial_refinement_level; // initial level of refinement
    unsigned int linear_solver_option; // linear solver option
    unsigned int scheme_option; // option for scheme used
    std::string entropy_string; // string for entropy function
    std::string entropy_derivative_string; // string for entropy derivative function
    double entropy_residual_coefficient; // value of entropy residual coefficient
    double jump_coefficient; // value of jump coefficient 
    // temporal discretization for entropy residual
    typename EntropyViscosity<dim>::TemporalDiscretization EV_time_discretization;
    bool output_meshes; // option to output meshes as .eps files
    bool is_steady_state; // is the problem steady-state?
    double end_time; // end time if transient problem is run
    double time_step_size; // time step size if transient problem is run
    bool output_exact_solution; // option to output exact solution
    bool output_initial_solution; // option to output initial solution
    bool output_DMP_bounds; // option to output DMP bounds
    double CFL_limit; // upper bound for the CFL number
    bool do_not_limit; // flag to turn off limiting for high order solution
    unsigned int n_quadrature_points; // number of quadrature points to use in formula
    bool save_convergence_results; // option to save convergence results
};

#include "TransportParameters.cc"
#endif
