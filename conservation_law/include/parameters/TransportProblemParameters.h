/**
 * \file TransportProblemParameters.h
 * \brief Provides the header for the TransportProblemParameters class.
 */
#ifndef TransportProblemParameters_h
#define TransportProblemParameters_h

#include <deal.II/base/parameter_handler.h>
#include "include/parameters/ProblemParameters.h"

using namespace dealii;

/**
 * \class TransportProblemParameters
 * \brief Class for parameters related to transport problems.
 */
template <int dim>
class TransportProblemParameters : public ProblemParameters<dim>
{
public:
  TransportProblemParameters(const std::string & problem_name);

  /** \brief Transport speed \f$v\f$ */
  double transport_speed;

  /** \brief Transport direction \f$\mathbf{\Omega}\f$ */
  Tensor<1, dim> transport_direction;

  /** \brief Function parser of cross section \f$\sigma\f$ */
  FunctionParser<dim> cross_section_function;

  /** \brief Function parser of source \f$q\f$ */
  FunctionParser<dim> source_function;

  /** \brief Flag that source is time-dependent */
  bool source_is_time_dependent;

protected:
  /** \brief Specification type for transport direction */
  std::string transport_direction_specification;

  /** \brief Azimuthal angle \f$\theta\f$ for transport direction */
  double azimuthal_angle;

  /** \brief Polar angle \f$\varphi\f$ for transport direction */
  double polar_angle;

  /** \brief x-component of transport direction, \f$\Omega_x\f$ */
  double transport_direction_x;

  /** \brief y-component of transport direction, \f$\Omega_y\f$ */
  double transport_direction_y;

  /** \brief z-component of transport direction, \f$\Omega_z\f$ */
  double transport_direction_z;

  /** \brief Option to normalize transport direction vector */
  bool normalize_transport_direction;

  /** \brief String for function parser of cross section \f$\sigma\f$ */
  std::string cross_section_string;

  /** \brief String for function parser of source \f$q\f$ */
  std::string source_string;

  /** \brief String for function parser of Dirichlet BC */
  std::string dirichlet_function_angularflux;

  /** \brief String for function parser of initial conditions */
  std::string initial_condition_angularflux;

  /** \brief String for function parser of exact solution */
  std::string exact_solution_angularflux;

  /** \brief Incoming flux value constant for function parsers */
  double incoming_value;

  /** \brief Cross section value 1 */
  double sigma1;
  /** \brief Cross section value 2 */
  double sigma2;

  /** \brief Source value constant for function parsers */
  double source_value;

  /** \brief x value 1 */
  double x1;
  /** \brief x value 2 */
  double x2;
  /** \brief x value 3 */
  double x3;
  /** \brief x value 4 */
  double x4;
  /** \brief y value 1 */
  double y1;
  /** \brief y value 2 */
  double y2;
  /** \brief y value 3 */
  double y3;

private:
  void declare_derived_parameters() override;

  void get_derived_parameters() override;

  void process_derived_parameters(
    Triangulation<dim> & triangulation,
    const FESystem<dim> & fe,
    const QGauss<dim - 1> & face_quadrature) override;

  void set_boundary_ids_incoming(Triangulation<dim> & triangulation,
                                 const FESystem<dim> & fe,
                                 const QGauss<dim - 1> & face_quadrature);
};

#include "src/parameters/TransportProblemParameters.cc"

#endif
