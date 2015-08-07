#ifndef EntropyViscosity_cc
#define EntropyViscosity_cc

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include "Viscosity.h"
#include "LowOrderViscosity.h"
#include "TransportParameters.h"

using namespace dealii;

/**
 * Class for computing entropy viscosity for a transport problem.
 */
template<int dim>
class EntropyViscosity : public Viscosity<dim>
{
public:

  EntropyViscosity(
    const FESystem<dim> &fe,
    const unsigned int &n_cells,
    const DoFHandler<dim> &dof_handler,
    const ConstraintMatrix &constraints,
    const QGauss<dim> &cell_quadrature,
    const QGauss<dim - 1> &face_quadrature,
    const Tensor<1, dim> &transport_direction,
    const FunctionParser<dim> &cross_section_function,
    FunctionParser<dim> &source_function,
    const std::string &entropy_string,
    const std::string &entropy_derivative_string,
    const double &entropy_residual_coefficient,
    const double &jump_coefficient,
    const double &domain_volume,
    const typename TransportParameters<dim>::TemporalDiscretization temporal_discretization,
    const LowOrderViscosity<dim> &low_order_viscosity,
    const SparseMatrix<double> &inviscid_matrix,
    SparseMatrix<double> &diffusion_matrix,
    SparseMatrix<double> &total_matrix);

  ~EntropyViscosity();

  void recomputeHighOrderSteadyStateMatrix(
    const Vector<double> & solution);

  void recompute_high_order_ss_matrix(
    const Vector<double> &old_solution,
    const Vector<double> &older_solution,
    const Vector<double> &oldest_solution,
    const double &old_dt,
    const double &older_dt,
    const double &time);

private:

  void compute_entropy_viscosity(
    const Vector<double> &old_solution,
    const Vector<double> &older_solution,
    const Vector<double> &oldest_solution,
    const double &old_dt,
    const double &older_dt,
    const double &time);
  void compute_normalization_constant(const Vector<double> &old_solution);
  void compute_temporal_discretization_constants(
    const double old_dt,
    const double older_dt);

  // mesh and dof data
  const FESystem<dim> *fe;
  const FEValuesExtractors::Scalar flux;
  const unsigned int n_dofs;
  const unsigned int faces_per_cell;

  // quadrature data
  const QGauss<dim> cell_quadrature;
  const QGauss<dim - 1> face_quadrature;
  const unsigned int n_q_points_cell;
  const unsigned int n_q_points_face;

  // physics data
  Tensor<1, dim> transport_direction;
  const FunctionParser<dim> *cross_section_function;
  FunctionParser<dim> *source_function;

  // entropy viscosity functions and data
  FunctionParser<1> entropy_function;
  FunctionParser<1> entropy_derivative_function;
  std::string entropy_string;
  std::string entropy_derivative_string;
  double entropy_residual_coefficient;
  double jump_coefficient;
  double domain_volume;
  double domain_averaged_entropy;
  double normalization_constant;

  // viscosity vectors
  Vector<double> entropy_viscosity;

  // temporal discretization for entropy residual
  typename TransportParameters<dim>::TemporalDiscretization temporal_discretization;

  // coefficients for entropy residual
  double a_old;
  double a_older;
  double a_oldest;
  double b_old;
  double b_older;

  // low-order viscosity
  const LowOrderViscosity<dim> * const low_order_viscosity;

  // matrices
  const SparseMatrix<double> * const inviscid_matrix;
  SparseMatrix<double> * const diffusion_matrix;
  SparseMatrix<double> * const total_matrix;
};

#include "EntropyViscosity.cc"
#endif
