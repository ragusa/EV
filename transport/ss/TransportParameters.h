#ifndef TransportParameters_cc
#define TransportParameters_cc

#include <deal.II/base/parameter_handler.h>

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
    bool use_adaptive_refinement; // option to use adaptive mesh refinement
    unsigned int n_refinement_cycles; // number of refinement cycles
    unsigned int initial_refinement_level; // initial level of refinement
    unsigned int linear_solver_option; // linear solver option
    bool output_solution; // option to output solution
    unsigned int n_quadrature_points; // number of quadrature points to use in formula
};

#include "TransportParameters.cc"
#endif
