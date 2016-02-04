#ifndef LowOrderViscosity_cc
#define LowOrderViscosity_cc

#include "Viscosity.h"
#include "PostProcessor.h"

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/constraint_matrix.h>

using namespace dealii;

/** \brief Class for low-order viscosity.
 */
template <int dim>
class LowOrderViscosity : public Viscosity<dim>
{
public:
  LowOrderViscosity(const unsigned int n_cells,
                    const unsigned int dofs_per_cell,
                    const DoFHandler<dim> & dof_handler,
                    const ConstraintMatrix & constraints,
                    const SparseMatrix<double> & inviscid_matrix,
                    SparseMatrix<double> & diffusion_matrix,
                    SparseMatrix<double> & total_matrix);
  ~LowOrderViscosity();
  /*
        void compute_bounds(const Vector<double>       &old_solution,
                            const SparseMatrix<double> &low_order_ss_matrix,
                            const Vector<double>       &ss_rhs,
                            const double               &dt);
  void output_bounds(const PostProcessor<dim> &postprocessor) const;
  */

private:
  void compute_viscous_bilinear_forms();
  void compute_low_order_viscosity(const SparseMatrix<double> & inviscid_matrix);

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> viscous_bilinear_forms;

  /*
        const DoFHandler<dim> *dof_handler;

        Vector<double> solution_min;
        Vector<double> solution_max;
  */
};

#include "LowOrderViscosity.cc"
#endif
