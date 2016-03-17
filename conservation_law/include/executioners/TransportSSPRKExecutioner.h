#ifndef TransportSSPRKExecutioner_cc
#define TransportSSPRKExecutioner_cc

#include "include/executioners/TransportTransientExecutioner.h"
#include "include/fct/TransportExplicitEulerFCT.h"
#include "include/parameters/RunParameters.h"
#include "include/time_integrators/SSPRKTimeIntegrator.h"

using namespace dealii;

/**
 * \brief Class for SSPRK executioner of a transport problem.
 */
template <int dim>
class TransportSSPRKExecutioner : public TransportTransientExecutioner<dim>
{
public:
  /** \brief Alias for scheme */
  using Scheme = typename TransportExecutioner<dim>::Scheme;
  /** \brief Alias for high-order scheme */
  using HighOrderScheme = typename TransportExecutioner<dim>::HighOrderScheme;

  TransportSSPRKExecutioner(const TransportRunParameters & parameters,
                            TransportProblemParameters<dim> & problem_parameters,
                            Triangulation<dim> & triangulation,
                            PostProcessor<dim> & postprocessor,
                            const double & nominal_dt);

protected:
  void compute_new_solution(const double & dt,
                            const double & dt_old,
                            const unsigned int & n) override;

private:
  void perform_fct_ssprk_step(
                              const double & dt,
                              const double & old_stage_dt,
                              const unsigned int & n);

  Vector<double> old_stage_solution;

  std::shared_ptr<TransportExplicitEulerFCT<dim>> fct;

  SSPRKTimeIntegrator<dim> ssprk;
};

#include "src/executioners/TransportSSPRKExecutioner.cc"
#endif
