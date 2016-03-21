/**
 * \file DMPThetaFCTFilter.h
 * \brief Provides the header for the DMPThetaFCTFilter class.
 */
#ifndef DMPThetaFCTFilter_h
#define DMPThetaFCTFilter_h

#include "include/fct/ThetaFCTFilter.h"

using namespace dealii;

/**
 * \brief Class for discrete maximum principle (DMP) filter for theta.
 */
template <int dim>
class DMPThetaFCTFilter : public ThetaFCTFilter<dim>
{
public:
  DMPThetaFCTFilter(const std::shared_ptr<Limiter<dim>> limiter,
                    const DoFHandler<dim> & dof_handler,
                    const FESystem<dim> & fe,
                    const SparseMatrix<double> & lumped_mass_matrix,
                    const double & theta);

protected:
  virtual void compute_solution_bounds(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs_new,
    const Vector<double> & ss_rhs_old) override;

  virtual void compute_antidiffusion_bounds(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs_new,
    const Vector<double> & ss_rhs_old,
    const Vector<double> & cumulative_antidiffusion) override;

  /** \brief temporary vector for minimum of new solution */
  Vector<double> solution_min_new;

  /** \brief temporary vector for maximum of new solution */
  Vector<double> solution_max_new;

  /** \brief temporary vector */
  Vector<double> tmp_vector;
};

#include "src/fct/DMPThetaFCTFilter.cc"

#endif
