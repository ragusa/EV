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

    void setup_system();
    void update_cell_sizes();
    void assemble_mass_matrix();
    void solve_erk();

    void update_flux_speeds();
    double compute_dt_from_cfl_condition();

    void mass_matrix_solve (Vector<double> &x);
    void linear_solve (const typename ConservationLawParameters<dim>::LinearSolverType     &linear_solver,
                       const SparseMatrix<double> &A,
                       const Vector<double>       &b,
                             Vector<double>       &x);

    void apply_Dirichlet_BC();

    // steady state residual functions
    void compute_ss_residual (Vector<double> &solution);
    virtual void compute_cell_ss_residual(FEValues<dim> &fe_values,
                                          const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          Vector<double> &cell_residual) = 0;
    virtual void compute_face_ss_residual(FEFaceValues<dim> &fe_face_values,
                                          const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          Vector<double> &cell_residual) = 0;
    virtual Tensor<1,dim> flux_derivative(const double u) = 0;
    void update_viscosities();
    void update_first_order_viscosities();
    void update_entropy_viscosities();
    void update_entropy_residuals();
    void update_jumps();
    virtual double entropy           (const double u) const = 0;
    virtual double entropy_derivative(const double u) const = 0;

    void output_results () const;

    void check_nan();

    virtual std::vector<std::string> get_component_names() = 0;
    virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
       get_component_interpretations() = 0;

    /** input parameters for conservation law */
    ConservationLawParameters<dim> conservation_law_parameters;
    /** number of components in the system */
    int n_components;

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

    /** quadrature formula for cells */
    const QGauss<dim>    quadrature;
    /** number of quadrature points for cells */
    const unsigned int   n_q_points_cell;
    /** quadrature formula for faces */
    const QGauss<dim-1>  face_quadrature;
    /** number of quadrature points for faces */
    const unsigned int   n_q_points_face;

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

    ConditionalOStream   verbose_cout;

    /** function parser for the initial condition expression given by user
     *  in input file */ 
    FunctionParser<dim>  initial_conditions;

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
    std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > entropy_residual_cell_q;
    std::map<typename DoFHandler<dim>::active_cell_iterator, Vector<double> > entropy_cell_q;
    std::map<typename DoFHandler<dim>::active_cell_iterator, double>          max_jumps_cell;
};

#include "ConservationLaw.cc"

#endif
