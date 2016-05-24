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
  DMPExplicitEulerFCTFilter(
    const RunParameters & run_parameters,
    const std::shared_ptr<Limiter<dim>> limiter,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const SparseMatrix<double> & lumped_mass_matrix_,
    const std::map<unsigned int, double> & dirichlet_values);

protected:
  virtual void compute_solution_bounds(const Vector<double> & old_solution,
                                       const double & dt,
                                       const Vector<double> & ss_reaction,
                                       const Vector<double> & ss_rhs,
                                       const double & t_old) override;
};

#include "src/fct/DMPExplicitEulerFCTFilter.cc"

#endif
