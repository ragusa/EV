/**
 * \file SteadyStateFCT.h
 * \brief Provides the header for the SteadyStateFCT class.
 */
#ifndef SteadyStateFCT_h
#define SteadyStateFCT_h

#include "include/fct/DMPSteadyStateFCTFilter.h"
#include "include/fct/FCT.h"

using namespace dealii;

/**
 * \brief Class for steady-state FCT.
 */
template <int dim>
class SteadyStateFCT : public FCT<dim>
{
public:
  SteadyStateFCT(const RunParameters & run_parameters,
                 const DoFHandler<dim> & dof_handler,
                 const FESystem<dim> & fe);

  void compute_antidiffusion_matrix(
    const Vector<double> & high_order_solution,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix);

  void compute_antidiffusion_vector(
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs,
    Vector<double> & antidiffusion_vector);

protected:
  void create_filters();

  virtual std::shared_ptr<SteadyStateFCTFilter<dim>> create_filter(
    const std::string & filter_string);

  void filter_antidiffusive_fluxes(
    const Vector<double> & solution,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix);

  /** \brief vector of FCT filters */
  std::vector<std::shared_ptr<SteadyStateFCTFilter<dim>>> filters;

  /** \brief flag to use cumulative antidiffusion algorithm */
  const bool use_cumulative_antidiffusion_algorithm;
};

#include "src/fct/SteadyStateFCT.cc"

#endif
