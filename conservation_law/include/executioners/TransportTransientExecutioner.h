#ifndef TransportTransientExecutioner_cc
#define TransportTransientExecutioner_cc

#include "include/executioners/TransportExecutioner.h"
#include "include/fct/TransportExplicitEulerFCT.h"
#include "include/parameters/RunParameters.h"
#include "include/time_integrators/SSPRKTimeIntegrator.h"

using namespace dealii;

/**
 * \brief Class for transient executioner of a transport problem.
 */
template <int dim>
class TransportTransientExecutioner : public TransportExecutioner<dim>
{
public:
  /** \brief Alias for scheme */
  using Scheme = typename TransportExecutioner<dim>::Scheme;
  /** \brief Alias for high-order scheme */
  using HighOrderScheme = typename TransportExecutioner<dim>::HighOrderScheme;
  /** \brief Alias for temporal discretization classification */
  using TemporalDiscretizationClassification =
    typename RunParameters::TemporalDiscretizationClassification;
  /** \brief Alias for time step size option */
  using TimeStepSizeOption = typename RunParameters::TimeStepSizeOption;

  TransportTransientExecutioner(
    const TransportRunParameters & parameters,
    TransportProblemParameters<dim> & problem_parameters,
    Triangulation<dim> & triangulation,
    PostProcessor<dim> & postprocessor,
    const double & nominal_dt);

  void run() override;

protected:
  virtual void compute_new_solution(const double & dt,
                                    const double & dt_old,
                                    const double & t_old,
                                    const unsigned int & n) = 0;

  SparseMatrix<double> consistent_mass_matrix;

  SparseMatrix<double> lumped_mass_matrix;

  Vector<double> old_solution;
  Vector<double> older_solution;
  Vector<double> oldest_solution;
  Vector<double> ss_rhs_new;

  const bool source_is_time_dependent;

  Vector<double> tmp_vector;

private:
  void assembleMassMatrices();

  double enforceCFLCondition(const double & dt_proposed) const;

  const TemporalDiscretizationClassification temporal_discretization;

  FunctionParser<dim> * const initial_conditions_function;

  double dt_nominal;

  double minimum_cell_diameter;
};

#include "src/executioners/TransportTransientExecutioner.cc"
#endif
