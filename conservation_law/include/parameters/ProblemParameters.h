/**
 * \file ProblemParameters.h
 * \brief Provides the header for the ProblemParameters class.
 */
#ifndef ProblemParameters_h
#define ProblemParameters_h

#include <string>
#include <sys/stat.h>

#include <deal.II/base/parameter_handler.h>

#include "include/bc/DirichletBoundaryConditions.h"
#include "include/other/CMakeVars.h"

using namespace dealii;

/**
 * \class ProblemParameters
 * \brief Class for conservation law problem parameters
 */
template <int dim>
class ProblemParameters
{
public:
  ProblemParameters(const std::string & problem_name,
                    const unsigned int & n_components,
                    const bool & specified_steady_state);

  void get_and_process_parameters(const std::string & parameters_file,
                                  Triangulation<dim> & triangulation,
                                  const FESystem<dim> & fe,
                                  const QGauss<dim - 1> & face_quadrature);

  /** \brief Problem name */
  const std::string problem_name;

  /** \brief Number of solution components */
  const unsigned int n_components;

  /** \brief Flag that problem is specified to be run in steady-state */
  const bool specified_steady_state;

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

  /** \brief Dirichlet BC function */
  std::shared_ptr<Function<dim>> dirichlet_function;

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
  /** \brief Flag that problem is transient */
  bool is_transient_problem;

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

  /** \brief Flag that exact solution should be used for Dirichlet BC */
  bool use_exact_solution_as_dirichlet_bc;

  /** \brief Flag that exact solution should be used for initial conditions */
  bool use_exact_solution_as_initial_conditions;

  /** \brief vector of initial condition function strings for each component */
  std::vector<std::string> initial_conditions_strings;

  /** \brief Type of function for exact solution */
  std::string exact_solution_type;

  /** \brief vector of exact solution function strings for each component */
  std::vector<std::string> exact_solution_strings;

  /** \brief vector of Dirichlet BC function strings for each component */
  std::vector<std::string> dirichlet_function_strings;

private:
  void declare_base_parameters();

  void get_base_parameters();

  void process_base_parameters();

  void process_shared_base_parameters(Triangulation<dim> & triangulation,
                                      const FESystem<dim> & fe,
                                      const QGauss<dim - 1> & face_quadrature);

  void generate_mesh_and_compute_volume(Triangulation<dim> & triangulation);

  void set_boundary_ids(Triangulation<dim> & triangulation,
                        const FESystem<dim> & fe,
                        const QGauss<dim - 1> & face_quadrature);

  virtual void declare_derived_parameters() = 0;

  virtual void get_derived_parameters() = 0;

  virtual void process_derived_parameters(
    Triangulation<dim> & triangulation,
    const FESystem<dim> & fe,
    const QGauss<dim - 1> & face_quadrature) = 0;
};

#include "src/parameters/ProblemParameters.cc"

#endif
