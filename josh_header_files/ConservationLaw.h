/** \file ConservationLaw.h
 *  \brief Provides the header for the ConservationLaw class.
 */
#ifndef ConservationLaw_h
#define ConservationLaw_h

#include <iostream>
#include <fstream>
#include <algorithm>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include "ConservationLawParameters.h"

using namespace dealii;

/** \class ConservationLaw
 *  \brief Class providing framework for solving a general conservation law.
 *
 *  This class solves a conservation law of the form:
 *  \f[
 *    \frac{\partial\vec{u}}{\partial t} 
 *    + \nabla \cdot \vec{f}(\vec{u}) = \vec{g}(\vec{u}),
 *  \f]
 *  where \f$\vec{u}\f$ is the vector of conservation variables, \f$\vec{f}(\vec{u})\f$
 *  is the flux, and \f$\vec{g}(\vec{u})\f$ is the forcing term.
 */
template <int dim>
class ConservationLaw
{
  public:

    ConservationLaw (const ConservationLawParameters<dim> &params);
    void run();

  protected:

    void initialize_system();
    void setup_system();
    void update_cell_sizes();
    void assemble_mass_matrix();
    void solve_runge_kutta();

    void update_flux_speeds();
    double compute_dt_from_cfl_condition();

    void mass_matrix_solve (Vector<double> &x);
    void linear_solve (const typename ConservationLawParameters<dim>::LinearSolverType     &linear_solver,
                       const SparseMatrix<double> &A,
                       const Vector<double>       &b,
                             Vector<double>       &x);

    void apply_Dirichlet_BC(const double &time);

    // steady state residual functions
    void compute_ss_residual (Vector<double> &solution);
    virtual void compute_cell_ss_residual(FEValues<dim> &fe_values,
                                          const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          Vector<double> &cell_residual) = 0;
    virtual void compute_face_ss_residual(FEFaceValues<dim> &fe_face_values,
                                          const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          Vector<double> &cell_residual) = 0;
    virtual Tensor<1,dim> flux_derivative(const double u) = 0;
    void update_viscosities(const double &dt);
    void update_first_order_viscosities();
    void update_entropy_viscosities(const double &dt);
    void update_entropy_residuals(const double &dt);
    void update_jumps();
    virtual double entropy           (const double u) const = 0;
    virtual double entropy_derivative(const double u) const = 0;

    void output_solution () const;
    void output_map(std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > &map,
                    const std::string &output_filename_base);
    void output_map(std::map<typename DoFHandler<dim>::active_cell_iterator, double> &map,
                    const std::string &output_filename_base);

    void check_nan();

    virtual std::vector<std::string> get_component_names() = 0;
    virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
       get_component_interpretations() = 0;

    virtual void define_problem() = 0;

    /** input parameters for conservation law */
    const ConservationLawParameters<dim> conservation_law_parameters;
    /** number of components in the system */
    const unsigned int n_components;

    /** triangulation; mesh */
    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;

    /** finite element system */
    const FESystem<dim>  fe;
    /** DoFs per cell */
    const unsigned int   dofs_per_cell;
    /** faces per cell */
    const unsigned int   faces_per_cell;
    /** DoF handler */
    DoFHandler<dim>      dof_handler;
    /** constraint matrix */
    ConstraintMatrix     constraints;

    /** number of quadrature points in each dimension */
    const unsigned int   n_q_points_per_dim;
    /** number of quadrature points per cell */
    const unsigned int   n_q_points_cell;
    /** number of quadrature points per face */
    const unsigned int   n_q_points_face;
    /** quadrature formula for cells */
    const QGauss<dim>    cell_quadrature;
    /** quadrature formula for faces */
    const QGauss<dim-1>  face_quadrature;

    /** solution of current time step */
    Vector<double>       current_solution;
    /** solution of previous time step */
    Vector<double>       old_solution;
    /** current time */
    double current_time;
    /** old time (beginning of each time step) */
    double old_time;
    /** system right-hand side */
    Vector<double>       system_rhs;

    /* sparsity pattern for the mass matrix */
    SparsityPattern      sparsity_pattern;
    /* mass matrix */
    SparseMatrix<double> mass_matrix;

    /** number of boundaries */
    unsigned int n_boundaries;
    /** enumeration for types of boundary conditions */
    enum BoundaryType {dirichlet};
    /** vector of types of boundary condition for each boundary indicator and component */
    std::vector<std::vector<BoundaryType> > boundary_types;
    /** vector of Dirichlet BC function strings, which will be parsed */
    std::vector<std::string>  dirichlet_function_strings;
    /** vector of Dirichlet BC functions created from parsed strings */
    FunctionParser<dim>       dirichlet_function;
    /** option to use exact solution function as Dirichlet BC */
    bool use_exact_solution_as_BC;
    /** option to skip computing the face residual if Dirichlet BCs used at
     *  all boundaries */
    bool need_to_compute_face_residual;

    /** initial conditions function strings for each component, which will be parsed */
    std::vector<std::string> initial_conditions_strings;
    /** initial conditions functions */
    FunctionParser<dim>      initial_conditions_function;

    /** option if the problem has an exact solution provided */
    bool has_exact_solution;
    /** exact solution function strings for each component, which will be parsed */
    std::vector<std::string> exact_solution_strings;
    /** exact solution functions */
    FunctionParser<dim>      exact_solution_function;

    /** vector of component names */
    std::vector<std::string> component_names;
    /** vector of data component interpretations (scalar or vector) */
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
       component_interpretations;

    /** minimum cell diameter; used for CFL condition */
    double dx_min;
    /** maximum flux speed; used for CFL condition */
    double max_flux_speed;
    /** max entropy deviation */
    double max_entropy_deviation;
    /** domain volume; used for calculation of domain-averaged entropy */
    double domain_volume;

    // maps
    std::map<typename DoFHandler<dim>::active_cell_iterator, double>          dx;
    std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > flux_speed_cell_q;
    std::map<typename DoFHandler<dim>::active_cell_iterator, double>          max_flux_speed_cell;
    std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > viscosity_cell_q;
    std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > first_order_viscosity_cell_q;
    std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > entropy_viscosity_cell_q;
    std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > entropy_viscosity_with_jumps_cell_q;
    std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > entropy_residual_cell_q;
    std::map<typename DoFHandler<dim>::active_cell_iterator, double>          max_entropy_residual_cell;
    std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > entropy_cell_q;
    std::map<typename DoFHandler<dim>::active_cell_iterator, double>          max_jumps_cell;
};

#include "ConservationLaw.cc"

#endif
