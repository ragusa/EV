/**
 * \file DMPExplicitEulerFCTFilter.h
 * \brief Provides the header for the DMPExplicitEulerFCTFilter class.
 */
#ifndef DMPExplicitEulerFCTFilter_h
#define DMPExplicitEulerFCTFilter_h

#include "include/fct/ExplicitEulerFCTFilter.h"

using namespace dealii;

/**
 * \brief Class for discrete maximum principle (DMP) filter for explicit Euler.
 */
template <int dim>
class DMPExplicitEulerFCTFilter : public ExplicitEulerFCTFilter<dim>
{
public:
  DMPExplicitEulerFCTFilter(const std::shared_ptr<Limiter> limiter,
                            const DoFHandler<dim> & dof_handler,
                            const SparseMatrix<double> & lumped_mass_matrix_);

protected:
  virtual void compute_solution_bounds(const Vector<double> & old_solution,
                                       const double & dt,
                                       const Vector<double> & ss_reaction,
                                       const Vector<double> & ss_rhs) override;

  virtual void compute_antidiffusion_bounds(
    const Vector<double> & old_solution,
    const double & dt,
    const Vector<double> & inviscid_ss_flux,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const Vector<double> & ss_rhs) override;

  /** \brief temporary vector */
  Vector<double> tmp_vector;
};

#include "src/fct/DMPExplicitEulerFCTFilter.cc"

#endif
