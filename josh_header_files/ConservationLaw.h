/*
 ConservationLaw class:
 modified from deal.II step-33 tutorial program
 
 solves a general conservation law of the form:
 
 du/dt + div(flux(c)) = rhs(c)
*/

#ifndef ConservationLaw_h
#define ConservationLaw_h

#include <iostream>
#include <fstream>

#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
//#include <deal.II/lac/sparse_direct.h>
//#include <deal.II/lac/solver_bicgstab.h>
//#include <deal.II/lac/precondition.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/vector_tools.h>
//#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include "ConservationLawParameters.h"

using namespace dealii;

template <int dim>
class ConservationLaw
{
  public:

    ConservationLaw (ParameterHandler &prm, const int &n_comp);//const std::string &input_file, const int &n_comp);
    void run ();

  private:

    void setup_system ();

    void assemble_system ();
    void assemble_cell_term (const FEValues<dim>             &fe_v,
                             const std::vector<unsigned int> &dofs);
    void assemble_face_term (const unsigned int               face_no,
                             const FEFaceValuesBase<dim>     &fe_v,
                             const FEFaceValuesBase<dim>     &fe_v_neighbor,
                             const std::vector<unsigned int> &dofs,
                             const std::vector<unsigned int> &dofs_neighbor,
                             const bool                       external_face,
                             const unsigned int               boundary_id,
                             const double                     face_diameter);

    std::pair<unsigned int, double> solve (Vector<double> &newton_update);

    void compute_refinement_indicators (Vector<double> &indicator) const;
    void refine_grid (const Vector<double> &indicator);

    void output_results () const;

    ConservationLawParameters<dim> conservation_law_parameters;
    int n_components;

    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;

    const FESystem<dim>  fe;
    DoFHandler<dim>      dof_handler;

    const QGauss<dim>    quadrature;
    const QGauss<dim-1>  face_quadrature;
	
    Vector<double>       current_solution;
    Vector<double>       old_solution;
    Vector<double>       old_old_solution;
    Vector<double>       predictor;

    Vector<double>       right_hand_side;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    ConditionalOStream   verbose_cout;

  protected:

    FunctionParser<dim>  initial_conditions;

    std::vector<std::string> component_names;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
       component_interpretations;
};

#include "ConservationLaw.cc"

#endif
