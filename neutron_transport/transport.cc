/* Author: Joshua Hansel
 * MATH 676 project, adapted from examples/step-9/step-9.cc
 * Description: solves the 1-d multi-group neutron transport equation
 */

/* preprocessor variable for number of dimensions:
 * will not compile for dim != 1 if this is not set to
 * be consistent with Hansel::dimension; exact solutions
 * are available only in 1-D
 */
#define DIMENSION 1

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/base/tensor_function.h>

#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/table.h>

#include <fstream>
#include <iostream>
#include <sstream>   // for stringstreams - creating file names
#include <cstdlib>   // for exit()
#include <algorithm> // for min()

namespace Hansel {

using namespace dealii;

const unsigned int dimension = DIMENSION; // number of spatial dimensions

/**
 * \class TransportProblem
 * \brief for defining a transport problem
 */
template<int dim>
class TransportProblem {
   public:
      // parameters class for defining input parameters
      class Parameters {
         public:
            Parameters();
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

      TransportProblem(const Parameters &parameters);
      ~TransportProblem();
      void run();

   private:
      void compute_directions();
      void setup_system();
      void assemble_system();
      void run_single();
      void solve();
      void refine_grid();
      void output_grid(const unsigned int cycle) const;
      void evaluate_error(const unsigned int cycle);
      void compute_viscous_bilinear_forms();
      void compute_max_principle_viscosity();
      void check_assembly_max_principle();
      void check_solution_max_principle();

      const Parameters &parameters; // input parameters
      unsigned int degree;
      unsigned int n_energy_groups;
      unsigned int n_directions;

      unsigned int n_variables; // total number of solution variables
      std::vector<Point<dim> > transport_directions; // direction vectors in Cartesian coordinates

      Triangulation<dim> triangulation;
      DoFHandler<dim> dof_handler;

      FESystem<dim> fe;
      const unsigned int dofs_per_cell;

      QGauss<dim>   cell_quadrature_formula;
      QGauss<dim-1> face_quadrature_formula;
      const unsigned int n_q_points_cell;
      const unsigned int n_q_points_face;

      ConstraintMatrix constraints;

      SparsityPattern constrained_sparsity_pattern;
      SparseMatrix<double> system_matrix;

      Vector<double> old_solution;
      Vector<double> present_solution;
      Vector<double> system_rhs;

      Vector<double> max_viscosity;
      Vector<double> entropy_viscosity;
      Vector<double> max_principle_viscosity;

      SparsityPattern unconstrained_sparsity_pattern;
      SparseMatrix<double> max_principle_viscosity_numerators;
      SparseMatrix<double> viscous_bilinear_forms;

      unsigned int nonlinear_iteration;

      ConvergenceTable convergence_table;

};

/**
 * \class ExactSolution1
 * \brief defines the exact solution for test problem 1
 */
template <int dim>
class ExactSolution1 : public Function<dim> {
   public:
      ExactSolution1 () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int component = 0) const;
};

/**
 * \fn    ExactSolution<1>::value
 * \brief computes the exact solution for test problem 1 for a given x-point
 */
template <>
double ExactSolution1<1>::value (const Point<1>  &p,
                                 const unsigned int) const
{
  const double sigma  = 1.0;
  const double source = 1.0 / (4.0 * numbers::PI);

  double return_value;
  return_value = source*(1.0 - std::exp(-sigma*(p[0]+1.0)));

  return return_value;
}

/**
 * \class ExactSolution2
 * \brief defines the exact solution for test problem 2
 */
template <int dim>
class ExactSolution2 : public Function<dim> {
   public:
      ExactSolution2 () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int component = 0) const;
};

/**
 * \fn    ExactSolution2<1>::value
 * \brief computes the exact solution for test problem 2 for a given x-point
 */
template <>
double ExactSolution2<1>::value (const Point<1>  &p,
                                 const unsigned int) const
{
  const double sigma = 100.0;
  const double incoming = 1.0;

  double return_value;
  if (p[0] < 0.0) return_value = incoming;
  else            return_value = incoming * std::exp(-sigma*p[0]);

  return return_value;
}

/**
 * \class ExactSolution4
 * \brief defines the exact solution for test problem 4
 */
template <int dim>
class ExactSolution4 : public Function<dim> {
   public:
      ExactSolution4 () : Function<dim>() {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int component = 0) const;
};

/**
 * \fn    ExactSolution4<1>::value
 * \brief computes the exact solution for test problem 4 for a given x-point
 */
template <>
double ExactSolution4<1>::value (const Point<1>  &p,
                                 const unsigned int) const
{
  double return_value = p[0] * std::sin(p[0]) * std::exp(-p[0]);
  return return_value;
}

/**
 * \class TotalSource
 * \brief defines the total source function
 */
template<int dim>
class TotalSource: public Function<dim> {
   public:
      TotalSource(const typename TransportProblem<dim>::Parameters &parameters) :
            Function<dim>(),
            parameters(parameters) {
      }

      double value(const unsigned int group,
            const unsigned int direction, const Point<dim> &p) const;

      void value_list(const unsigned int group,
            const unsigned int direction,
            const std::vector<Point<dim> > &points,
            std::vector<double> &values) const;
   private:
      const typename TransportProblem<dim>::Parameters &parameters;
};

/**
 * \fn    TotalSource<dim>::value
 * \brief computes the total source at a given point in space
 */
template<int dim>
double TotalSource<dim>::value(const unsigned int group,
      const unsigned int direction, const Point<dim> &p) const {
   Assert(group < 2, ExcNotImplemented());
   double return_value = 0.0;
   switch (parameters.source_option) {
      case 1: {
         if (group == 0)
            return_value = parameters.source_value / (4.0 * numbers::PI); // isotropic source term
         else
            return_value = 0.0 / (4.0 * numbers::PI); // isotropic source term
         break;
      } case 2: {
         if (group == 0) {
            bool in_nonzero_region = true;
            for (unsigned int i = 0; i < dimension; ++i)
               if (p[i] < 0.0) {
                  in_nonzero_region = false;
                  break;
               }
            if (in_nonzero_region)
               return_value = parameters.source_value;
         } else
            return_value = 0.0;
         break;
      } case 3: {
         /* manufactured solution source
          * solution is: x*sin(x)*exp(-x) and sigma = 1
          */
         return_value = std::sin(p[0]) * std::exp(-p[0])
                         + p[0] * std::cos(p[0]) * std::exp(-p[0]);
         break;
      }
      default: {
         Assert(false,ExcNotImplemented());
         break;
      }
   }
   return return_value;
}

/**
 * \fn    TotalSource<dim>::value_list
 * \brief computes the total source at a number of points in space
 */
template<int dim>
void TotalSource<dim>::value_list(const unsigned int group,
      const unsigned int direction, const std::vector<Point<dim> > &points,
      std::vector<double> &values) const {
   Assert(values.size() == points.size(),
         ExcDimensionMismatch (values.size(), points.size()));

   for (unsigned int i = 0; i < points.size(); ++i)
      values[i] = value(group, direction, points[i]);
}

/**
 * \class TotalCrossSection
 * \brief defines the total cross section function
 */
template<int dim>
class TotalCrossSection: public Function<dim> {
   public:
      TotalCrossSection(const typename TransportProblem<dim>::Parameters &parameters) :
            Function<dim>(),
            parameters(parameters) {
      }

      double value(const unsigned int group,
            const unsigned int direction, const Point<dim> &p) const;

      void value_list(const unsigned int group,
            const unsigned int direction,
            const std::vector<Point<dim> > &points,
            std::vector<double> &values) const;
   private:
      const typename TransportProblem<dim>::Parameters parameters;
};

/**
 * \fn    TotalCrossSection<dim>::value
 * \brief computes the total cross section at a given point in space
 */
template<int dim>
double TotalCrossSection<dim>::value(const unsigned int group,
      const unsigned int direction, const Point<dim> &p) const {
   Assert(group < 2, ExcNotImplemented());
   double return_value = 0.0;
   switch (parameters.total_cross_section_option) {
      case 1: {
         if (group == 0)
            return_value = parameters.total_cross_section_value;
         else if (group == 1)
            return_value = 0.0;
         else {
            Assert(false,ExcNotImplemented());
         }
         break;
      }
      case 2: {
         if (group == 0) {
            bool in_nonzero_region = true;
            for (unsigned int i = 0; i < dimension; ++i)
               if (p[i] < 0.0) {
                  in_nonzero_region = false;
                  break;
               }
            if (in_nonzero_region)
               return_value = parameters.total_cross_section_value;
         } else if (group == 1)
            return_value = 0.0;
         else {
            Assert(false,ExcNotImplemented());
         }
         break;
      }
      default: {
         Assert(false,ExcNotImplemented());
         break;
      }
   }

   return return_value;
}

/**
 * \fn    TotalCrossSection<dim>::value_list
 * \brief computes the total cross section at a number of points in space
 */
template<int dim>
void TotalCrossSection<dim>::value_list(const unsigned int group,
      const unsigned int direction, const std::vector<Point<dim> > &points,
      std::vector<double> &values) const {
   Assert(values.size() == points.size(),
         ExcDimensionMismatch (values.size(), points.size()));

   for (unsigned int i = 0; i < points.size(); ++i)
      values[i] = TotalCrossSection<dim>::value(group, direction, points[i]);
}

/**
 * \fn    TransportProblem<dim>::TransportProblem
 * \brief constructor for TransportProblem class
 */
template<int dim>
TransportProblem<dim>::TransportProblem(const Parameters &parameters) :
      parameters(parameters),
      degree(parameters.degree),
      n_energy_groups(parameters.n_energy_groups),
      n_directions(parameters.n_directions),
      n_variables(n_energy_groups * n_directions),
      dof_handler(triangulation),
      fe(FE_Q<dim>(degree), n_variables),
      dofs_per_cell(fe.dofs_per_cell),
      cell_quadrature_formula(degree+1),
      face_quadrature_formula(degree+1),
      n_q_points_cell(cell_quadrature_formula.size()),
      n_q_points_face(face_quadrature_formula.size()),
      nonlinear_iteration(0)
{
   // compute the direction unit vectors
   compute_directions();
}

/**
 * \fn    TransportProblem<dim>::~TransportProblem
 * \brief destructor for TransportProblem class
 */
template<int dim>
TransportProblem<dim>::~TransportProblem() {
   dof_handler.clear();
}

/**
 * \fn    TransportProblem<dim>::Parameters::Parameters
 * \brief constructor for the TransportProblem<dim>::Parameters class
 */
template<int dim>
TransportProblem<dim>::Parameters::Parameters() :
      degree(1), n_energy_groups(1), n_directions(1), use_adaptive_mesh_refinement(
            false), n_refinement_cycles(5), initial_refinement_level(2), solver_option(
            1), preconditioner_option(1), source_option(1), source_value(0.0e0), total_cross_section_option(
            1), total_cross_section_value(1.0e0), incoming_flux(1.0e0), viscosity_type(0),
            max_viscosity_coefficient(5.0e-1), entropy_viscosity_coefficient(1.0e-1),
            max_nonlinear_iterations(10), relative_difference_tolerance(1.0e-6),
            exact_solution_id(0), output_meshes(false)
            {
}

/**
 * \fn    TransportProblem<dim>::Parameters::declare_parameters
 * \brief defines all of the input parameters
 */
template<int dim>
void TransportProblem<dim>::Parameters::declare_parameters(
      ParameterHandler &prm) {
   prm.declare_entry("Finite element degree", "1", Patterns::Integer(),
         "Polynomial degree of finite elements");
   prm.declare_entry("Number of energy groups", "1", Patterns::Integer(),
         "Number of neutron energy groups");
   prm.declare_entry("Number of directions", "1", Patterns::Integer(),
         "Number of transport directions");
   prm.declare_entry("Use adaptive mesh refinement", "false",
         Patterns::Bool(),
         "Option to use adaptive mesh refinement instead of uniform");
   prm.declare_entry("Number of refinement cycles", "5", Patterns::Integer(),
         "Number of mesh refinement cycles");
   prm.declare_entry("Initial refinement level", "2", Patterns::Integer(),
         "Number of refinements for first mesh refinement cycle");
   prm.declare_entry("Solver option", "1", Patterns::Integer(),
         "Option for linear solver");
   prm.declare_entry("Preconditioner option", "1", Patterns::Integer(),
         "Option for preconditioner for linear solver");
   prm.declare_entry("Source option", "1", Patterns::Integer(),
         "Option for source definition");
   prm.declare_entry("Source value", "0.0e0", Patterns::Double(),
         "Value of extraneous source term");
   prm.declare_entry("Total cross section option", "1", Patterns::Integer(),
         "Option for total cross section definition");
   prm.declare_entry("Total cross section value", "1.0e0", Patterns::Double(),
         "Value of total cross section");
   prm.declare_entry("Incoming flux", "1.0e0", Patterns::Double(),
         "Value for flux on incoming boundary");
   prm.declare_entry("Viscosity type", "0", Patterns::Integer(),
         "Option for viscosity type: none, first-order, or entropy");
   prm.declare_entry("Max viscosity coefficient", "5.0e-1", Patterns::Double(),
         "Coefficient for the first-order viscosity");
   prm.declare_entry("Entropy viscosity coefficient", "1.0e-1", Patterns::Double(),
         "Coefficient for the entropy viscosity");
   prm.declare_entry("Maximum number of nonlinear iterations", "10", Patterns::Integer(),
         "Maximum number of nonlinear iterations to use when using entropy viscosity");
   prm.declare_entry("Relative difference tolerance", "1.0e-6", Patterns::Double(),
         "Relative difference tolerance for nonlinear convergence");
   prm.declare_entry("Exact solution ID", "0", Patterns::Integer(),
         "ID of exact solution to use when evaluating error");
   prm.declare_entry("Output mesh", "false", Patterns::Bool(),
         "Option to output meshes as .eps files");
}

/**
 * \fn    TransportProblem<dim>::Parameters::get_parameters
 * \brief get the input parameters
 */
template<int dim>
void TransportProblem<dim>::Parameters::get_parameters(ParameterHandler &prm) {
   degree = prm.get_integer("Finite element degree");
   n_energy_groups = prm.get_integer("Number of energy groups");
   n_directions = prm.get_integer("Number of directions");
   use_adaptive_mesh_refinement = prm.get_bool("Use adaptive mesh refinement");
   n_refinement_cycles = prm.get_integer("Number of refinement cycles");
   initial_refinement_level = prm.get_integer("Initial refinement level");
   solver_option = prm.get_integer("Solver option");
   preconditioner_option = prm.get_integer("Preconditioner option");
   source_option = prm.get_integer("Source option");
   source_value = prm.get_double("Source value");
   total_cross_section_option = prm.get_integer("Total cross section option");
   total_cross_section_value = prm.get_double("Total cross section value");
   incoming_flux = prm.get_double("Incoming flux");
   viscosity_type = prm.get_integer("Viscosity type");
   max_viscosity_coefficient = prm.get_double(
         "Max viscosity coefficient");
   entropy_viscosity_coefficient = prm.get_double(
         "Entropy viscosity coefficient");
   max_nonlinear_iterations = prm.get_integer("Maximum number of nonlinear iterations");
   relative_difference_tolerance = prm.get_double("Relative difference tolerance");
   exact_solution_id = prm.get_integer("Exact solution ID");
   output_meshes = prm.get_bool("Output mesh");
}

/**
 * \fn    TransportProblem<dim>::compute_directions
 * \brief creates the transport directions
 */
template<int dim>
void TransportProblem<dim>::compute_directions() {
   Assert(dim < 3, ExcNotImplemented());
   // radians between each transport direction
   double angle_interval = 2 * numbers::PI / n_directions;
   // first angle is always 0 radians
   double angle_radians = 0.0;
   for (unsigned int k = 0; k < n_directions; ++k) {
      // unit vector of direction in Cartesian coordinates - restricted
      //  to 2 dimensions
      Point<dim> direction(true);
      if (dim == 1)
         direction(0) = std::cos(angle_radians);
      else if (dim == 2) {
         direction(0) = std::cos(angle_radians);
         direction(1) = std::sin(angle_radians);
      } else
         Assert(false,ExcNotImplemented());

      // put direction unit vector in vector of all directions
      transport_directions.push_back(direction);
      // increment angle
      angle_radians += angle_interval;
   }
}

/**
 * \fn    TransportProblem<dim>::setup_system
 * \brief set up the problem before assembly of the linear system
 */
template<int dim>
void TransportProblem<dim>::setup_system()
{
   dof_handler.distribute_dofs(fe);

   // reinitialize viscosity vectors
   max_viscosity          .reinit(triangulation.n_active_cells());
   entropy_viscosity      .reinit(triangulation.n_active_cells());
   max_principle_viscosity.reinit(triangulation.n_active_cells());

   // clear constraint matrix and make hanging node constraints for new mesh
   constraints.clear();
   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
   constraints.close();

   // reinitialize sparsity pattern of system matrix
   CompressedSparsityPattern compressed_constrained_sparsity_pattern(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern(dof_handler,
                                               compressed_constrained_sparsity_pattern,
                                               constraints,
                                               false);
   constrained_sparsity_pattern.copy_from(compressed_constrained_sparsity_pattern);

   // allocate matrices to be used with maximum-principle preserving viscosity
   // reinitialize sparsity pattern of auxiliary matrices
   CompressedSparsityPattern compressed_unconstrained_sparsity_pattern(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern(dof_handler, compressed_unconstrained_sparsity_pattern);
   unconstrained_sparsity_pattern.copy_from(compressed_unconstrained_sparsity_pattern);

   // reinitialize auxiliary matrices with sparsity pattern
   viscous_bilinear_forms            .reinit(unconstrained_sparsity_pattern);
   max_principle_viscosity_numerators.reinit(unconstrained_sparsity_pattern);

   // compute viscous bilinear forms
   compute_viscous_bilinear_forms();

   // reinitialize solution vector, system matrix, and rhs
   system_matrix   .reinit(constrained_sparsity_pattern);
   present_solution.reinit(dof_handler.n_dofs());
   system_rhs      .reinit(dof_handler.n_dofs());
}

/** \fn void TransportProblem<dim>::compute_viscous_bilinear_forms()
 *  \brief Computes viscous bilinear forms, to be used in the computation of
 *         maximum-principle preserving first order viscosity.
 *
 *         Each element of the resulting matrix, \f$B_{i,j}\f$ is computed as
 *         follows:
 *         \f[
 *            B_{i,j} = \sum_{K\subset S_{i,j}}b_K(\varphi_i,\varphi_j)
 *         \f]
 */
template <int dim>
void TransportProblem<dim>::compute_viscous_bilinear_forms()
{
   viscous_bilinear_forms = 0; // zero out matrix

   // dofs per cell
   unsigned int dofs_per_cell = fe.dofs_per_cell;

   // local dof indices
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // query cell volume
      double cell_volume = cell->measure();

      // add local bilinear forms to global matrix
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
         for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            double b_cell = 0.0;
            if (j == i) {
               b_cell = cell_volume;
            } else {
               b_cell = -1.0/(dofs_per_cell-1.0)*cell_volume;
            }
            viscous_bilinear_forms.add(local_dof_indices[i],
                                       local_dof_indices[j],
                                       b_cell);
         }
      }
   }
}

/**
 * \fn    TransportProblem<dim>::assemble_system
 * \brief assemble the system matrix and right hand side
 */
template<int dim>
void TransportProblem<dim>::assemble_system() {
   max_principle_viscosity_numerators = 0;

   const TotalCrossSection<dim> total_cross_section(parameters);
   const TotalSource<dim> total_source(parameters);

   // FE values, for assembly terms
   FEValues<dim> fe_values(fe, cell_quadrature_formula,
         update_values | update_gradients | update_quadrature_points
               | update_JxW_values);

   // FE face values, for boundary conditions
   FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
         update_values | update_quadrature_points | update_JxW_values
               | update_normal_vectors | update_gradients);

   FullMatrix<double> inviscid_cell_matrix(dofs_per_cell, dofs_per_cell);
   FullMatrix<double>          cell_matrix(dofs_per_cell, dofs_per_cell);
   Vector<double> cell_rhs(dofs_per_cell);

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // total cross section values for an energy group and direction at each
   //  quadrature point on cell
   std::vector<double> total_cross_section_values(n_q_points_cell);
   // total source values for an energy group and direction at each
   //  quadrature point on cell
   std::vector<double> total_source_values(n_q_points_cell);

   // create extractors
   Table<2, FEValuesExtractors::Scalar*> fluxes(n_energy_groups, n_directions);
   for (unsigned int g = 0; g < n_energy_groups; ++g)
      for (unsigned int d = 0; d < n_directions; ++d)
         fluxes(g, d) = new FEValuesExtractors::Scalar(g * n_directions + d);

   // cell iterator
   typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active(), endc = dof_handler.end();

   // compute domain volume for denominator of domain-averaged entropy
   double domain_volume = 0.0;
   for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
      // reinitialize FE values
      fe_values.reinit(cell);
      // loop over quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         domain_volume += fe_values.JxW(q);
   }

   // loop over energy groups
   for (unsigned int g = 0; g < n_energy_groups; ++g) {
      // loop over transport directions
      for (unsigned int d = 0; d < n_directions; ++d) {
         // get domain-averaged entropy if using entropy viscosity for this iteration
         // ---------------------------------------------------------------------
         double domain_averaged_entropy;
         double max_entropy_deviation_domain = 0.0;
         if ((parameters.viscosity_type == 2)&&(nonlinear_iteration != 0)) {
            // compute domain-averaged entropy
            double domain_integral_entropy = 0.0;
            // loop over cells
            for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
               // reinitialize FE values
               fe_values.reinit(cell);
               // loop over quadrature points
               for (unsigned int q = 0; q < n_q_points_cell; ++q) {
                  // get angular flux at quadrature point
                  double angular_flux = 0.0;
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                     angular_flux += fe_values[*(fluxes(g, d))].value(j, q);
                  }
                  // compute entropy at quadrature point
                  double entropy = 1.0/2.0 * angular_flux * angular_flux;
                  // add contribution of quadrature point to entropy integral
                  domain_integral_entropy += entropy * fe_values.JxW(q);
               }
            }
            // domain-averaged entropy
            domain_averaged_entropy = domain_integral_entropy / domain_volume;

            // find max deviation of entropy from domain entropy average
            // loop over cells
            for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
               // reinitialize FE values
               fe_values.reinit(cell);
               // get old values and gradients
               std::vector<double> old_values(n_q_points_cell);
               fe_values[*(fluxes(g, d))].get_function_values(old_solution,old_values);
               // loop over quadrature points
               for (unsigned int q = 0; q < n_q_points_cell; ++q) {
                  // get angular flux at quadrature point
                  double angular_flux = 0.0;
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                     angular_flux += old_values[q];
                  }
                  // compute entropy at quadrature point
                  double entropy = 1.0/2.0 * angular_flux * angular_flux;
                  // add contribution of quadrature point to entropy integral
                  max_entropy_deviation_domain = std::max(max_entropy_deviation_domain,
                        std::abs(entropy-domain_averaged_entropy));
               }
            }
         }

         // loop over cells
         unsigned int i_cell = 0;
         for (cell = dof_handler.begin_active(); cell != endc; ++cell, ++i_cell) {
            // initialize local matrix and rhs to zero
            inviscid_cell_matrix = 0;
            cell_matrix = 0;
            cell_rhs = 0;

            // get local dof indices
            cell->get_dof_indices(local_dof_indices);

            // reinitialize FE values
            fe_values.reinit(cell);

            // get total cross section for all quadrature points
            total_cross_section.value_list(g, d,
                  fe_values.get_quadrature_points(),
                  total_cross_section_values);
            // get total source for all quadrature points
            total_source.value_list(g, d, fe_values.get_quadrature_points(),
                  total_source_values);

            // compute viscosity
            // ------------------------------------------------------------------
            double viscosity;
            double h = cell->diameter();
            double max_viscosity_cell = parameters.max_viscosity_coefficient * h;
            max_viscosity(i_cell) = max_viscosity_cell;
            switch (parameters.viscosity_type) {
               case 0: {
                  break;
               } case 1: {
                  viscosity = max_viscosity_cell;
                  break;
               } case 2: {
                  if (nonlinear_iteration != 0) {
                     // get old values and gradients
                     std::vector<double> old_values(n_q_points_cell);
                     std::vector<Tensor<1,dim> > old_gradients(n_q_points_cell);
                     fe_values[*(fluxes(g, d))].get_function_values(old_solution,old_values);
                     fe_values[*(fluxes(g, d))].get_function_gradients(old_solution,old_gradients);

                     // compute entropy values at each quadrature point on cell
                     std::vector<double> entropy_values(n_q_points_cell,0.0);
                     for (unsigned int q = 0; q < n_q_points_cell; ++q)
                        entropy_values[q] = 1.0/2.0 * old_values[q] * old_values[q];
                     // compute entropy residual values at each quadrature point on cell
                     std::vector<double> entropy_residual_values(n_q_points_cell,0.0);
                     for (unsigned int q = 0; q < n_q_points_cell; ++q)
                        entropy_residual_values[q] = std::abs(transport_directions[d] *
                              old_values[q] * old_gradients[q]
                                                              + total_cross_section_values[q] * entropy_values[q]);
                     // compute entropy deviation values at each quadrature point on cell
                     std::vector<double> entropy_deviation_values(n_q_points_cell,0.0);
                     for (unsigned int q = 0; q < n_q_points_cell; ++q)
                        entropy_deviation_values[q] = std::abs(entropy_values[q] - domain_averaged_entropy);
                     // determine cell maximum entropy residual and entropy deviation
                     double max_entropy_residual = 0.0;
                     for (unsigned int q = 0; q < n_q_points_cell; ++q) {
                        max_entropy_residual = std::max(max_entropy_residual,entropy_residual_values[q]);
                     }

                     // compute entropy viscosity
                     double entropy_viscosity_cell = parameters.entropy_viscosity_coefficient *
                           h * h * max_entropy_residual / max_entropy_deviation_domain;
                     entropy_viscosity(i_cell) = entropy_viscosity_cell;

                     // determine viscosity: minimum of first-order viscosity and entropy viscosity
                     viscosity = std::min(max_viscosity_cell, entropy_viscosity_cell);
                  } else
                     viscosity = max_viscosity_cell;

                  break;
               } case 3: {
                  // maximum-principle preserving viscosity does not add Laplacian term;
                  // to avoid unnecessary branching in assembly loop, just add Laplacian
                  // term with zero viscosity, so just keep viscosity with its initialization
                  // of zero
                  viscosity = 0.0;
                  break;
               } default: {
                  Assert(false,ExcNotImplemented());
                  break;
               }
            }
            // compute cell contributions to global system
            // ------------------------------------------------------------------
            for (unsigned int q = 0; q < n_q_points_cell; ++q) {
               for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                     // store integrals of divergence and total interaction term
                     // so that they may be used in computation of max-principle
                     // preserving viscosity
                     inviscid_cell_matrix(i,j) += (
                        // divergence term
                        fe_values[*(fluxes(g, d))].value(i, q)
                           * transport_directions[d]
                           * fe_values[*(fluxes(g, d))].gradient(j, q) +
                        // total interaction term
                        fe_values[*(fluxes(g, d))].value(i, q)
                           * total_cross_section_values[q]
                           * fe_values[*(fluxes(g, d))].value(j, q)
                        ) * fe_values.JxW(q);

                     // add to matrix
                     cell_matrix(i, j) +=
                        // viscosity term
                        viscosity
                           * fe_values[*(fluxes(g, d))].gradient(i, q)
                           * fe_values[*(fluxes(g, d))].gradient(j, q)
                           * fe_values.JxW(q);
                  } // end j

                  cell_rhs(i) +=
                     // total source term
                     fe_values[*(fluxes(g, d))].value(i, q)
                        * total_source_values[q] * fe_values.JxW(q);
               } // end i
            } // end q

            // add face terms; these arise due to integration by parts of viscosity term
            // ------------------------------------------------------------------
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
            {
               if (cell->face(face)->at_boundary()) {
                  fe_face_values.reinit(cell, face);

                  for (unsigned int q = 0; q < n_q_points_face; ++q)
                     for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        for (unsigned int j = 0; j < dofs_per_cell; ++j)
                           cell_matrix(i, j) -= (viscosity
                                 * fe_face_values.shape_value(i, q)
                                 * fe_face_values.normal_vector(q)
                                 * fe_face_values.shape_grad(j, q)
                                 * fe_face_values.JxW(q));
               }
            } 

            // add inviscid terms to cell matrix
            cell_matrix.add(1.0, inviscid_cell_matrix);

            // aggregate local matrix and rhs to global matrix and rhs
            constraints.distribute_local_to_global(cell_matrix,
                                                   cell_rhs,
                                                   local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);

            // aggregate local inviscid matrix global matrix
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
               for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  max_principle_viscosity_numerators.add(local_dof_indices[i],
                                                         local_dof_indices[j],
                                                         inviscid_cell_matrix(i,j));
                  
         } // end cell
      } // end d
   } // end g

   // add viscous bilinear form for maximum-principle preserving viscosity
   // ---------------------------------------------------------------------------
   if (parameters.viscosity_type == 3) {
      // compute maximum-principle preserving viscosity
      compute_max_principle_viscosity();

      for (unsigned int g = 0; g < n_energy_groups; ++g) {
         for (unsigned int d = 0; d < n_directions; ++d) {
            unsigned int i_cell = 0;
            for (cell = dof_handler.begin_active(); cell != endc; ++cell, ++i_cell) {
               // reset cell matrix to zero
               cell_matrix = 0;

               // compute cell volume
               double cell_volume = cell->measure();
               // compute cell contribution to global system matrix
               for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                     double viscous_bilinear_form;
                     if (j == i)
                        viscous_bilinear_form = cell_volume;
                     else
                        viscous_bilinear_form = -1.0/(dofs_per_cell - 1.0)*cell_volume;
 
                     cell_matrix(i,j) += max_principle_viscosity(i_cell) * viscous_bilinear_form;
                  }
               }

               // aggregate local matrix and rhs to global matrix and rhs
               cell->get_dof_indices(local_dof_indices);
               constraints.distribute_local_to_global(cell_matrix,
                                                      local_dof_indices,
                                                      system_matrix);

            }
         }
      }
   }

   // apply boundary conditions
   // ---------------------------------------------------------------------------
   // loop over transport directions
   for (unsigned int d = 0; d < n_directions; ++d) {
      // reset boundary indicators to zero
      for (cell = dof_handler.begin_active(); cell != endc; ++cell)
         for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               ++face)
            if (cell->face(face)->at_boundary())
               cell->face(face)->set_boundary_indicator(0);

      // loop over cells
      for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
         // loop over faces of cell
         for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               ++face) {
            // if face is at boundary
            if (cell->face(face)->at_boundary()) {
               // reinitialize FE face values
               fe_face_values.reinit(cell, face);
               // determine if the transport flux is incoming through this face;
               //  it isn't necessary to loop over all face quadrature points because
               //  the transport direction and normal vector are the same at each
               //  quadrature point; therefore, quadrature point 0 is arbitrarily chosen
               double small = -1.0e-12;
               if (fe_face_values.normal_vector(0) * transport_directions[d] < small) {
                  // mark boundary as incoming flux boundary: indicator 1
                  cell->face(face)->set_boundary_indicator(1);
               }
            }
         }
      }
      // create component masks for current transport direction
      std::vector<bool> component_mask(n_variables, false);
      // loop over energy groups
      for (unsigned int g = 0; g < n_energy_groups; ++g)
         component_mask[g * n_directions + d] = true;
      // apply Dirichlet boundary condition
      std::map<unsigned int, double> boundary_values;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               ConstantFunction<dim>(parameters.incoming_flux, n_variables),
                                               boundary_values,
                                               component_mask);
      MatrixTools::apply_boundary_values(boundary_values,
                                         system_matrix,
                                         present_solution,
                                         system_rhs);
   } // end transport directions

   // deallocate extractors
   for (unsigned int g = 0; g < n_energy_groups; ++g)
      for (unsigned int d = 0; d < n_directions; ++d)
         delete fluxes(g, d);

} // end assembly

/** \fn void TransportProblem<dim>::compute_max_principle_viscosity()
 *  \brief Computes the maximum-principle preserving first order viscosity for each cell.
 */
template <int dim>
void TransportProblem<dim>::compute_max_principle_viscosity()
{
   // local dof indices
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   // loop over cells to compute first order viscosity at each quadrature point
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   unsigned int i_cell = 0;
   for (; cell != endc; ++cell, ++i_cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         max_principle_viscosity(i_cell) = 0.0;
         for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
               if (i != j) {
                  max_principle_viscosity(i_cell) = std::max(max_principle_viscosity(i_cell),
                     std::abs(max_principle_viscosity_numerators(local_dof_indices[i],local_dof_indices[j]))/
                     (-viscous_bilinear_forms(local_dof_indices[i],local_dof_indices[j])));
               }
            }
         }
      }
   }
}

/**
 * \fn    TransportProblem<dim>::solve
 * \brief solve the linear system
 */
template<int dim>
void TransportProblem<dim>::solve() {
   switch (parameters.solver_option) {
      case 1: {
         SparseDirectUMFPACK A_direct;
         A_direct.initialize(system_matrix);
         A_direct.vmult(present_solution, system_rhs);
         break;
      }
      case 2: {
         SolverControl solver_control(1000, 1e-6);
         SolverBicgstab<> solver(solver_control);

         switch (parameters.preconditioner_option) {
            case 1: {
               solver.solve(system_matrix, present_solution, system_rhs,
                     PreconditionIdentity());
               break;
            }
            case 2: {
               PreconditionJacobi<> preconditioner;
               preconditioner.initialize(system_matrix, 1.0);
               solver.solve(system_matrix, present_solution, system_rhs,
                     preconditioner);
               break;
            }
            default: {
               Assert(false,ExcNotImplemented());
               break;
            }
         }
         break;
      }
      default: {
         Assert(false,ExcNotImplemented());
         break;
      }
   }

   constraints.distribute(present_solution);
}

/**
 * \fn    TransportProblem<dim>::refine_grid
 * \brief refine the grid
 */
template<int dim>
void TransportProblem<dim>::refine_grid() {
   if (parameters.use_adaptive_mesh_refinement) {
      Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

      KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim - 1>(3),
            typename FunctionMap<dim>::type(), present_solution,
            estimated_error_per_cell);

      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
            estimated_error_per_cell, 0.3, 0.03);

      triangulation.execute_coarsening_and_refinement();
   } else
      triangulation.refine_global(1);
}

/**
 * \fn    TransportProblem<dim>::output_grid
 * \brief output the grid of the given cycle
 */
template<int dim>
void TransportProblem<dim>::output_grid(const unsigned int cycle) const {
   // create grid output file name
   std::string filename = "grid-";
   filename += ('0' + cycle);
   Assert(cycle < 10, ExcInternalError());

   filename += ".eps";
   std::ofstream output(filename.c_str());

   // write grid to eps
   GridOut grid_out;
   grid_out.write_eps(triangulation, output);
}

/**
 * \fn    TransportProblem<dim>::run
 * \brief run the problem
 */
template<int dim>
void TransportProblem<dim>::run() {
   // loop over refinement cycles
   for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; ++cycle) {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      // generate mesh if in first cycle, else refine
      if (cycle == 0) {
         GridGenerator::hyper_cube(triangulation, -1, 1);
         triangulation.refine_global(parameters.initial_refinement_level);
      } else {
         refine_grid();
      }

      // setup system - distribute finite elements, reintialize matrices and vectors
      setup_system();

      std::cout << "   Number of active cells:       "
            << triangulation.n_active_cells();
      std::cout      << std::endl;
      std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs();
      std::cout      << std::endl;

      // solve system. if using entropy viscosity, begin nonlinear iteration, else
      // just solve linear system
      nonlinear_iteration = 0;
      if (parameters.viscosity_type == 2) {
         bool converged = false; // converged nonlinear iteration
         for (unsigned int iter = 0; iter < parameters.max_nonlinear_iterations; ++iter) {
            std::cout << "   Nonlinear iteration " << iter;
            if (iter == 0) std::cout << std::endl;
            assemble_system();
            solve();
            // if not the first iteration, evaluate the convergence criteria
            if (nonlinear_iteration != 0) {
               // evaluate the difference between the current and previous solution iterate
               double old_norm = old_solution.l2_norm();
               old_solution -= present_solution;
               double difference_norm = old_solution.l2_norm();
               double relative_difference = difference_norm / old_norm;
               std::cout << ": Error: " << relative_difference << std::endl;
               if (relative_difference < parameters.relative_difference_tolerance) {
                  converged = true;
                  break;
               }
            }
            // update the old solution and iteration number
            old_solution = present_solution;
            nonlinear_iteration++;
         }
         // report if the solution did not converge
         if (!converged) {
            std::cout << "The solution did not converge in " << parameters.max_nonlinear_iterations << " iterations";
            std::cout << std::endl;
         }
      } else {
         // system is linear and requires just one solve
         assemble_system();
         //check_assembly_max_principle();
         solve();
         check_solution_max_principle();
      }

      // evaluate errors
      if (parameters.exact_solution_id != 0) evaluate_error(cycle);

      // output grid
      if (parameters.output_meshes) {
         if (dim > 1)
            output_grid(cycle);
      }
   }

   // output and post-processing
   //============================================================================

   // create the name string for each flux variable
   std::vector<std::string> variable_names;
   // loop over neutron energy groups
   for (unsigned int g = 0; g < n_energy_groups; ++g) {
      // loop over transport directions
      for (unsigned int k = 0; k < n_directions; ++k) {
         // variable name includes energy group and transport direction index
         std::stringstream variable_stringstream;
         variable_stringstream << "flux_group" << g << "_direction" << k;
         std::string variable_name = variable_stringstream.str();
         variable_names.push_back(variable_name);
      }
   }
   // create output data object
   DataOut<dim> data_out;
   data_out.attach_dof_handler(dof_handler);
   data_out.add_data_vector(present_solution, variable_names);
   data_out.build_patches(degree + 1);

   // create output filename
   std::string output_extension;
   if (dim == 1) output_extension = ".gpl";
   else          output_extension = ".vtk";
   std::string viscosity_string;
   switch (parameters.viscosity_type) {
      case 0: {
         viscosity_string = "none";
         break;
      } case 1: {
         viscosity_string = "first_order";
         break;
      } case 2: {
         viscosity_string = "entropy";
         break;
      } case 3: {
         viscosity_string = "max_principle";
         break;
      } default: {
         Assert(false, ExcNotImplemented());
         break;
      }
   }
   std::stringstream output_filename_ss;
   output_filename_ss << "output/solution_" << viscosity_string;
   if (parameters.exact_solution_id != 0) output_filename_ss << "_" << parameters.exact_solution_id;
   output_filename_ss << output_extension;
   std::string output_filename = output_filename_ss.str();
   char *output_filename_char = (char*)output_filename.c_str();
   std::ofstream output(output_filename_char);

   // if 1-d, write solution for gnuplot; otherwise, vtk
   if (dim == 1) data_out.write_gnuplot(output);
   else          data_out.write_vtk(output);

   // write viscosity plot files
   if (parameters.viscosity_type != 0) {
      DataOut<dim> visc_out;
      visc_out.attach_dof_handler(dof_handler);
      // add viscosity data vector(s)
      if ((parameters.viscosity_type == 1)||(parameters.viscosity_type == 2))
         visc_out.add_data_vector(max_viscosity,"Max_Viscosity",DataOut<dim>::type_cell_data);
      if (parameters.viscosity_type == 2)
         visc_out.add_data_vector(entropy_viscosity,"Entropy_Viscosity",DataOut<dim>::type_cell_data);
      if (parameters.viscosity_type == 3)
         visc_out.add_data_vector(max_principle_viscosity,"Max_Principle_Viscosity",DataOut<dim>::type_cell_data);
      // build patches and write to file
      visc_out.build_patches(degree + 1);
      if (dim == 1) {
         std::ofstream visc_out_stream("output/viscosity.gpl");
         visc_out.write_gnuplot(visc_out_stream);
      } else {
         std::ofstream visc_out_stream("output/viscosity.vtk");
         visc_out.write_vtk(visc_out_stream);
      }
   }

   // print convergence table
   if (parameters.exact_solution_id != 0) {
      convergence_table.set_precision("cell size", 3);
      convergence_table.set_scientific("cell size", true);
      convergence_table.set_precision("L2 error", 3);
      convergence_table.set_scientific("L2 error", true);
      convergence_table.evaluate_convergence_rates("L2 error", ConvergenceTable::reduction_rate_log2);
      convergence_table.write_text(std::cout);
   }
}

/**
 * \fn    TransportProblem<dim>::evaluate_error
 * \brief evaluate error between numerical and exact solution
 */
template<int dim>
void TransportProblem<dim>::evaluate_error(const unsigned int cycle) {
   Assert(dim == 1, ExcNotImplemented());

   // error per cell
   Vector<double> difference_per_cell (triangulation.n_active_cells());

   // compute error with analytic solution
   switch (parameters.exact_solution_id) {
      case 1: { // Test problem 1
         ExactSolution1<1> exact_solution;

#if DIMENSION == 1
         VectorTools::integrate_difference (MappingQ<1>(1), dof_handler, present_solution,
               exact_solution,
               difference_per_cell,
               QGauss<1>(degree+1),
               VectorTools::L2_norm);
#endif
         break;
      } case 2: { // Test problem 2
         ExactSolution2<1> exact_solution;

#if DIMENSION == 1
         VectorTools::integrate_difference (MappingQ<1>(1), dof_handler, present_solution,
               exact_solution,
               difference_per_cell,
               QGauss<1>(degree+1),
               VectorTools::L2_norm);
#endif
         break;
      } case 4: {
         ExactSolution4<1> exact_solution;

#if DIMENSION == 1
         VectorTools::integrate_difference (MappingQ<1>(1), dof_handler, present_solution,
               exact_solution,
               difference_per_cell,
               QGauss<1>(degree+1),
               VectorTools::L2_norm);
#endif
         break;
      } default: {
         Assert(false, ExcNotImplemented())
         break;
      }
   }

   // compute L2 error of vector of cell errors
   const double L2_error = difference_per_cell.l2_norm();

   const unsigned int n_active_cells = triangulation.n_active_cells();
   const unsigned int n_dofs = dof_handler.n_dofs();
   const double avg_cell_length = std::pow(2.0,dim) / std::pow(n_active_cells,1.0/dim);

   // add error values to convergence table
   convergence_table.add_value("cycle", cycle);
   convergence_table.add_value("cells", n_active_cells);
   convergence_table.add_value("dofs", n_dofs);
   convergence_table.add_value("cell size", avg_cell_length);
   convergence_table.add_value("L2 error", L2_error);
}

template<int dim>
void TransportProblem<dim>::check_assembly_max_principle()
{
   // check value of maximum principle preserving viscosity
   for (unsigned int i = 0; i < triangulation.n_active_cells(); ++i)
      std::cout << "nu[" << i << "] = " << max_principle_viscosity(i) << std::endl;

   // check system matrix
   std::ofstream matrix_out("output/matrix.txt");
   system_matrix.print_formatted(matrix_out,2,true,0,"0",1);
   matrix_out.close();

   // check denominators of viscosity
   std::ofstream denom_out("output/denom.txt");
   viscous_bilinear_forms.print_formatted(denom_out,2,true,0,"0",1);
   denom_out.close();

   // check numerators of viscosity
   std::ofstream num_out("output/num.txt");
   max_principle_viscosity_numerators.print_formatted(num_out,2,true,0,"0",1);
   num_out.close();
}

// check that solution at a node is bounded by its neighbors
template<int dim>
void TransportProblem<dim>::check_solution_max_principle()
{
   // local dof indices
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   // loop over cells to compute first order viscosity at each quadrature point
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   unsigned int i_cell = 0;
   for (; cell != endc; ++cell, ++i_cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
      }
   }
}

}

/**
 * \fn    main
 * \brief reads input file and then runs problem
 */
int main(int argc, char ** argv) {
   try {
      dealii::deallog.depth_console(0);

      // input filename
      std::string input_filename;
      if (argc < 2) {
         std::cout << "Need to provide an argument to specify the input file.\n"
            << "The argument '#' specifies the input file 'test_problem_#.in'"
            << std::endl;
         std::exit(1);
      } else {
         std::stringstream input_filename_ss;
         input_filename_ss << "input/test_problem_" << argv[1] << ".in";
         input_filename = input_filename_ss.str();
      }

      // get input parameters
      dealii::ParameterHandler parameter_handler;
      Hansel::TransportProblem<Hansel::dimension>::Parameters parameters;
      parameters.declare_parameters(parameter_handler);
      parameter_handler.read_input(input_filename);
      parameters.get_parameters(parameter_handler);

      // run problem
      Hansel::TransportProblem<Hansel::dimension> transport_problem(parameters);
      transport_problem.run();

   } catch (std::exception &exc) {
      std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what()
            << std::endl << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;
      return 1;
   } catch (...) {
      std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;
      return 1;
   };

   return 0;
}
