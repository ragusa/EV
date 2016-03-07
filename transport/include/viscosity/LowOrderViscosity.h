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

private:
  void compute_viscous_bilinear_forms();
  void compute_low_order_viscosity(const SparseMatrix<double> & inviscid_matrix);

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> viscous_bilinear_forms;
};

#include "LowOrderViscosity.cc"
#endif
