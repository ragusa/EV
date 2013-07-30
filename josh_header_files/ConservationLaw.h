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
//#include <deal.II/lac/precondition.h>
//#include <deal.II/lac/compressed_sparsity_pattern.h>

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

    ConservationLaw (ParameterHandler &prm, const int &n_comp);
    void run ();

  private:

    void solve_erk();
    void setup_system ();

    void linear_solve (const typename ConservationLawParameters<dim>::LinearSolverType     &linear_solver,
                       const SparseMatrix<double> &A,
                       const Vector<double>       &b,
                             Vector<double>       &x);
    void invert_mass_matrix (const Vector<double> &b, Vector<double> &x);
    virtual void compute_ss_residual (double t, Vector<double> &solution) = 0;

    void output_results () const;

    /** input parameters for conservation law */
    ConservationLawParameters<dim> conservation_law_parameters;
    int n_components;

    /** polynomial degree of finite elements */
    int degree;

    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;

    const FESystem<dim>  fe;
    DoFHandler<dim>      dof_handler;

    const QGauss<dim>    quadrature;
    const QGauss<dim-1>  face_quadrature;

    Vector<double>       current_solution;
    Vector<double>       old_solution;
    Vector<double>       ss_residual;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;

    ConditionalOStream   verbose_cout;

  protected:

    FunctionParser<dim>  initial_conditions;

    std::vector<std::string> component_names;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
       component_interpretations;
};

#include "ConservationLaw.cc"

#endif
