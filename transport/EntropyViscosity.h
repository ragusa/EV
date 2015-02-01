#ifndef EntropyViscosity_cc
#define EntropyViscosity_cc

#include <deal.II/base/quadrature_lib.h>
//#include <deal.II/base/logstream.h>
#include <deal.II/base/function_parser.h>
//#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/vector.h>
//#include <deal.II/lac/full_matrix.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
//#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/dofs/dof_handler.h>
//#include <deal.II/dofs/dof_accessor.h>
//#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
//#include <deal.II/fe/mapping_q.h>

//#include <deal.II/numerics/vector_tools.h>
//#include <deal.II/numerics/matrix_tools.h>

//#include <fstream>
//#include <iostream>
//#include <sstream>
//#include <cstdlib>
//#include <algorithm>

using namespace dealii;

/** \brief Class for computing entropy viscosity for a transport problem
 */
template<int dim>
class EntropyViscosity {
   public:
      EntropyViscosity(const FESystem<dim>   &fe,
                       const unsigned int    &n_cells,
                       const DoFHandler<dim> &dof_handler,
                       const QGauss<dim>     &cell_quadrature,
                       const QGauss<dim-1>   &face_quadrature,
                       const Tensor<1,dim>   &transport_direction,
                       const FunctionParser<dim> &cross_section_function,
                       const FunctionParser<dim> &source_function,
                       const std::string     &entropy_string,
                       const double          &entropy_residual_coefficient,
                       const double          &jump_coefficient,
                       const double          &domain_volume);
      ~EntropyViscosity();
      Vector<double> compute_entropy_viscosity(const Vector<double> &new_solution,
                                               const Vector<double> &old_solution,
                                               const double         &dt);

   private:
      void compute_entropy_domain_average(const Vector<double> &new_solution,
                                          const Vector<double> &old_solution);

      // mesh and dof data
      const FESystem<dim> *fe;
      const FEValuesExtractors::Scalar flux;
      const unsigned int n_cells;
      const DoFHandler<dim> *dof_handler;
      const unsigned int n_dofs;
      const unsigned int dofs_per_cell;
      const unsigned int faces_per_cell;

      // quadrature data
      const QGauss<dim>   cell_quadrature;
      const QGauss<dim-1> face_quadrature;
      const unsigned int n_q_points_cell;
      const unsigned int n_q_points_face;

      // physics data
      Tensor<1,dim> transport_direction;
      const FunctionParser<dim> *cross_section_function;
      const FunctionParser<dim> *source_function;

      // entropy viscosity functions and data
      FunctionParser<dim> entropy_function;
      std::string entropy_string;
      double entropy_residual_coefficient;
      double jump_coefficient;
      double domain_volume;
      double domain_averaged_entropy;
      double max_entropy_deviation_domain;

      // viscosity vectors
      Vector<double> entropy_viscosity;
};

#include "EntropyViscosity.cc"
#endif
