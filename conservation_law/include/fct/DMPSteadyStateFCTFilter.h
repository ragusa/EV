/**
 * \file DMPSteadyStateFCTFilter.h
 * \brief Provides the header for the DMPSteadyStateFCTFilter class.
 */
#ifndef DMPSteadyStateFCTFilter_h
#define DMPSteadyStateFCTFilter_h

#include "include/fct/SteadyStateFCTFilter.h"

using namespace dealii;

/**
 * \brief Class for discrete maximum principle (DMP) filter for steady-state.
 */
template <int dim>
class DMPSteadyStateFCTFilter : public SteadyStateFCTFilter<dim>
{
public:
  DMPSteadyStateFCTFilter(const RunParameters & run_parameters,
                          const std::shared_ptr<Limiter<dim>> limiter,
                          const DoFHandler<dim> & dof_handler,
                          const FESystem<dim> & fe);

protected:
  virtual void compute_solution_bounds(
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs) override;

  virtual void compute_antidiffusion_bounds(
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs,
    const Vector<double> & cumulative_antidiffusion) override;
};

#include "src/fct/DMPSteadyStateFCTFilter.cc"

#endif
