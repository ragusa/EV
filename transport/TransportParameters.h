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
    unsigned int degree; // polynomial degree of finite elements
    unsigned int n_energy_groups; // number of energy groups
    unsigned int n_directions; // number of discrete ordinates
    bool use_adaptive_mesh_refinement; // option to use adaptive mesh refinement
    unsigned int n_refinement_cycles; // number of refinement cycles
    unsigned int initial_refinement_level; // initial level of refinement
    unsigned int solver_option; // solver option
    unsigned int preconditioner_option; // preconditioner option
    unsigned int source_option; // source option
    double source_value; // maximum source value
    unsigned int total_cross_section_option; // total cross section option
    double total_cross_section_value; // maximum total cross section value
    double incoming_flux; // value for incoming flux
    unsigned int viscosity_type; // option for viscosity type
    double max_viscosity_coefficient; // value of maximum viscosity coefficient
    double entropy_viscosity_coefficient; // value of entropy viscosity coefficient
    unsigned int max_nonlinear_iterations; // maximum number of nonlinear iterations
    double relative_difference_tolerance; // relative difference tolerance for nonlinear convergence
    unsigned int exact_solution_id; // ID of exact solution to use when evaluating error
    bool output_meshes; // option to output meshes as .eps files
};

#include "TransportParameters.cc"
#endif
