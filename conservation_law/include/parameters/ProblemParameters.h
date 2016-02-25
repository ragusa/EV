/**
 * \file ProblemParameters.h
 * \brief Provides the header for the ProblemParameters class.
 */
#ifndef ProblemParameters_h
#define ProblemParameters_h

#include <string>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

/**
 * \class ProblemParameters
 * \brief Class for conservation law problem parameters
 */
template <int dim>
class ProblemParameters
{
public:
  ProblemParameters();

  void get_and_process_parameters(
  const std::string & filename, ParameterHandler & parameter_handler,
                                  Triangulation<dim> & triangulation);

  /** \brief initial conditions function */
  FunctionParser<dim> initial_conditions_function;

  /** \brief Flag that problem has an exact solution provided */
  bool has_exact_solution;

  /** \brief Exact solution function */
  std::shared_ptr<Function<dim>> exact_solution_function;

  /** \brief type of boundary conditions */
  std::string boundary_conditions_type;

  /** \brief boundary conditions */
  std::shared_ptr<BoundaryConditions<dim>> boundary_conditions;

  /** \brief vector of Dirichlet BC functions for each Dirichlet boundary */
  std::vector<std::shared_ptr<FunctionParser<dim>>> dirichlet_function;

  /** \brief Domain volume */
  double domain_volume;

  /** \brief Flag that test problem has a default end time */
  bool has_default_end_time;

  /** \brief Default end time for test problem */
  double default_end_time;

protected:

  /** \brief Parameter handler */
  ParameterHandler parameter_handler;

  /** \brief constants for function parsers */
  std::map<std::string, double> constants;

  /** \brief Flag that problem is valid in 1-D */
  bool valid_in_1d;
  /** \brief Flag that problem is valid in 2-D */
  bool valid_in_2d;
  /** \brief Flag that problem is valid in 3-D */
  bool valid_in_3d;

  /** \brief Shape description of domain */
  std::string domain_shape;

  /** \brief For box-type domains, the starting x position */
  double x_start;
  /** \brief For box-type domains, the starting y position */
  double y_start;
  /** \brief For box-type domains, the starting z position */
  double z_start;

  /** \brief For box-type domains, the domain width in x direction */
  double x_width;
  /** \brief For box-type domains, the domain width in y direction */
  double y_width;
  /** \brief For box-type domains, the domain width in z direction */
  double z_width;

  /** \brief Scheme for determination of boundary IDs for each face */
  std::string boundary_id_scheme;

  /** \brief Type of boundary conditions */
  std::string boundary_conditions_type;

  /** \brief Flag that exact solution should be used for Dirichlet BC */
  bool use_exact_solution_as_dirichlet_bc;

  /** \brief number of Dirichlet boundaries */
  unsigned int n_dirichlet_boundaries;

  /** \brief vector of initial condition function strings for each component */
  std::vector<std::string> initial_condition_strings;

  /** \brief Type of function for exact solution */
  std::string exact_solution_type;

  /** \brief vector of exact solution function strings for each component */
  std::vector<std::string> exact_solution_strings;

  /** \brief vector of Dirichlet BC function strings for each Dirichlet
             boundary and component */
  std::vector<std::vector<std::string>> dirichlet_function_strings;

private:

  void declare_base_parameters();

  void get_base_parameters();

  void process_base_parameters(Triangulation<dim> & triangulation);

  void generate_mesh_and_compute_volume(Triangulation<dim> & triangulation);

  void set_boundary_ids(Triangulation<dim> & triangulation);

  virtual void declare_derived_parameters() = 0;

  virtual void get_derived_parameters() = 0;

  virtual void process_derived_parameters() = 0;
};

#include "src/parameters/ProblemParameters.cc"

#endif
