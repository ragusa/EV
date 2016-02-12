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
  TransportProblemParameters();

  static void declare_parameters(ParameterHandler & parameter_handler);

  void get_parameters(ParameterHandler & parameter_handler);

  /** \brief Transport speed \f$v\f$ */
  double transport_speed;

  /** \brief x-component of transport direction, \f$\Omega_x\f$ */
  double transport_direction_x;

  /** \brief y-component of transport direction, \f$\Omega_y\f$ */
  double transport_direction_y;

  /** \brief z-component of transport direction, \f$\Omega_z\f$ */
  double transport_direction_z;

  /** \brief String for function parser of cross section \f$\sigma\f$ */
  std::string cross_section_string;

  /** \brief String for function parser of source \f$q\f$ */
  std::string source_string;

  /** \brief Incoming flux value constant for function parsers */
  double incoming_value;

  /** \brief Cross section value constant for function parsers */
  double cross_section_value;

  /** \brief Source value constant for function parsers */
  double source_value;

  std::string dirichlet_function;
  std::string initial_condition;
  std::string exact_solution;
};

#include "src/parameters/TransportProblemParameters.cc"

#endif
