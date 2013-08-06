/** \file ConservationLaw.h
 *  \brief Provides the header for the ConservationLaw class.
 */
#ifndef ConservationLaw_h
#define ConservationLaw_h

#include <iostream>
#include <fstream>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

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

    void solve_erk();
    void setup_system();
    void linear_solve (const typename ConservationLawParameters<dim>::LinearSolverType     &linear_solver,
                       const SparseMatrix<double> &A,
                       const Vector<double>       &b,
                             Vector<double>       &x);
    void invert_mass_matrix (const Vector<double> &b, Vector<double> &x);
    virtual void compute_ss_residual (double t, Vector<double> &solution) = 0;
    void output_results () const;
    virtual std::vector<std::string> get_component_names() = 0;
    virtual std::vector<DataComponentInterpretation::DataComponentInterpretation>
       get_component_interpretations() = 0;
    double compute_dt_from_cfl_condition();
    void check_nan();

    /** input parameters for conservation law */
    ConservationLawParameters<dim> conservation_law_parameters;
    /** number of components in the system */
    int n_components;

    /** triangulation; mesh */
    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;

    /** finite element system */
    const FESystem<dim>  fe;
    /** DoF handler */
    DoFHandler<dim>      dof_handler;

    /** quadrature formula for cells */
    const QGauss<dim>    quadrature;
    /** quadrature formula for faces */
    const QGauss<dim-1>  face_quadrature;

    /** solution of current time step */
    Vector<double>       current_solution;
    /** solution of previous time step */
    Vector<double>       old_solution;
    /** steady-state residual; used in intermediate stages */
    Vector<double>       ss_residual;

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
};

#include "ConservationLaw.cc"

#endif
