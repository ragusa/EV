#ifndef Viscosity_cc
#define Viscosity_cc

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>

using namespace dealii;

/**
 * Class for a general cell-based artificial viscosity.
 */
template <int dim>
class Viscosity
{
public:
  Viscosity(const unsigned int n_cells,
            const unsigned int dofs_per_cell,
            const DoFHandler<dim> & dof_handler,
            const ConstraintMatrix & constraints);
  ~Viscosity();

  void add_diffusion_matrix(const SparseMatrix<double> & inviscid_matrix,
                            SparseMatrix<double> & diffusion_matrix,
                            SparseMatrix<double> & total_matrix);
  double get_viscosity_value(const unsigned int i) const;
  void output_viscosity(const std::string output_file) const;

protected:
  void compute_diffusion_matrix(SparseMatrix<double> & diffusion_matrix);

  Vector<double> viscosity;

  const unsigned int n_cells;
  const unsigned int n_dofs;
  const unsigned int dofs_per_cell;
  const DoFHandler<dim> * const dof_handler;
  const ConstraintMatrix * const constraints;
};

#include "Viscosity.cc"
#endif
