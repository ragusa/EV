#ifndef EntropyViscosity_cc
#define EntropyViscosity_cc

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_parser.h>

#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

using namespace dealii;

/** \brief Class for computing entropy viscosity for a transport problem
 */
template<int dim>
class EntropyViscosity {
   public:
      // temporal discretization type for entropy residual
      enum TemporalDiscretization {FE, BE, CN, BDF2};

      EntropyViscosity(const FESystem<dim>   &fe,
                       const unsigned int    &n_cells,
                       const DoFHandler<dim> &dof_handler,
                       const QGauss<dim>     &cell_quadrature,
                       const QGauss<dim-1>   &face_quadrature,
                       const Tensor<1,dim>   &transport_direction,
                       const FunctionParser<dim> &cross_section_function,
                       FunctionParser<dim>   &source_function,
                       const std::string     &entropy_string,
                       const std::string     &entropy_derivative_string,
                       const double          &entropy_residual_coefficient,
                       const double          &jump_coefficient,
                       const double          &domain_volume,
                       const TemporalDiscretization temporal_discretization);

      ~EntropyViscosity();

      Vector<double> compute_entropy_viscosity(const Vector<double> &old_solution,
                                               const Vector<double> &older_solution,
                                               const Vector<double> &oldest_solution,
                                               const double         &old_dt,
                                               const double         &older_dt,
                                               const double         &time);

   private:
      void compute_normalization_constant(const Vector<double> &old_solution);
      void compute_temporal_discretization_constants(const double old_dt, const double older_dt);

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
      FunctionParser<dim> *source_function;

      // entropy viscosity functions and data
      FunctionParser<dim> entropy_function;
      FunctionParser<dim> entropy_derivative_function;
      std::string entropy_string;
      std::string entropy_derivative_string;
      double entropy_residual_coefficient;
      double jump_coefficient;
      double domain_volume;
      double domain_averaged_entropy;
      double normalization_constant;

      // viscosity vectors
      Vector<double> entropy_viscosity;

      TemporalDiscretization temporal_discretization;

      // coefficients for entropy residual
      double a_old;
      double a_older;
      double a_oldest;
      double b_old;
      double b_older;
};

#include "EntropyViscosity.cc"
#endif
