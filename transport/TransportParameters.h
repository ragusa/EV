#ifndef TransportParameters_cc
#define TransportParameters_cc

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/** \brief parameters class for defining input parameters
*/
//template <int dim>
class TransportParameters {
 public:
    TransportParameters();
    static void declare_parameters(ParameterHandler &prm);
    void get_parameters(ParameterHandler &prm);

    unsigned int problem_id; // problem ID
    unsigned int degree; // polynomial degree of finite elements
    unsigned int time_integrator_option; // time integrator option
    bool use_adaptive_mesh_refinement; // option to use adaptive mesh refinement
    unsigned int n_refinement_cycles; // number of refinement cycles
    unsigned int initial_refinement_level; // initial level of refinement
    unsigned int solver_option; // solver option
    unsigned int preconditioner_option; // preconditioner option
    unsigned int scheme_option; // option for scheme used
    double old_first_order_viscosity_coefficient; // value of old first order viscosity coefficient
    double entropy_viscosity_coefficient; // value of entropy viscosity coefficient
    double jump_coefficient; // value of jump coefficient
    bool output_meshes; // option to output meshes as .eps files
    bool is_steady_state; // is the problem steady-state?
    double end_time; // end time if transient problem is run
    double time_step_size; // time step size if transient problem is run
    bool output_exact_solution; // option to output exact solution
    bool output_initial_solution; // option to output initial solution
    bool output_DMP_bounds; // option to output DMP bounds
    double CFL_limit; // upper bound for the CFL number
    bool do_not_limit; // flag to turn off limiting for high order solution
};

#include "TransportParameters.cc"
#endif
