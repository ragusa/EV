/**
 * \file TransportExactDMPSSFCTFilter.h
 * \brief Provides the header for the TransportExactDMPSSFCTFilter class.
 */
#ifndef TransportExactDMPSSFCTFilter_h
#define TransportExactDMPSSFCTFilter_h

#include "include/fct/DMPSteadyStateFCTFilter.h"
#include "include/fct/SteadyStateFCTFilter.h"
#include "include/parameters/TransportProblemParameters.h"

using namespace dealii;

/**
 * \brief Class to evaluate the DMP solution bounds \f$W_i^-\f$ and \f$W_i^+\f$
 *        with the exact solution for steady-state.
 *
 * These solution bounds are not practical, as they require the exact solution,
 * but these bounds are useful as an academic exercise to determine if steady-
 * state FCT convergence difficulties are resolved when the solution bounds
 * are no longer implicit.
 *
 * The DMP solution bounds are computed as
 * \f[
 *   W_i^-(\mathrm{U}) = -\frac{1}{A^L_{i,j}}\sum\limits_{j\ne i}
 *     A^L_{i,j}U_{min,i} + \frac{b_i}{A^L_{i,i}} \,,
 * \f]
 * \f[
 *   W_i^+(\mathrm{U}) = -\frac{1}{A^L_{i,j}}\sum\limits_{j\ne i}
 *     A^L_{i,j}U_{max,i} + \frac{b_i}{A^L_{i,i}} \,,
 * \f]
 * where instead of
 * \f[
 *   U_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} U_j
 * \f]
 * \f[
 *   U_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} U_j
 * \f]
 * the exact solution is used in place of solution values:
 * \f[
 *   U_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} u^{exact}(\mathbf{x}_j)
 * \f]
 * \f[
 *   U_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} u^{exact}(\mathbf{x}_j)
 * \f]
 */
template <int dim>
class TransportExactDMPSSFCTFilter : public SteadyStateFCTFilter<dim>
{
public:
  TransportExactDMPSSFCTFilter(
    const RunParameters & run_parameters,
    const std::shared_ptr<Limiter<dim>> limiter,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const std::map<unsigned int, double> & dirichlet_values,
    TransportProblemParameters<dim> & problem_parameters);

  void compute_solution_bounds(const Vector<double> & solution,
                               const SparseMatrix<double> & low_order_ss_matrix,
                               const Vector<double> & ss_rhs) override;

protected:
  /** \brief DMP filter */
  DMPSteadyStateFCTFilter<dim> dmp_filter;

  /** \brief exact solution values at each node */
  Vector<double> exact_solution_values;
};

#include "src/fct/TransportExactDMPSSFCTFilter.cc"

#endif
