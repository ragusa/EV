#ifndef Viscosity_cc
#define Viscosity_cc

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

using namespace dealii;

/** \brief Class for a general cell-based artificial viscosity.
 */
template<int dim>
class Viscosity {
   public:
      Viscosity(const unsigned int      n_cells,
                const unsigned int      dofs_per_cell,
                const DoFHandler<dim>  &dof_handler);
      ~Viscosity();
      void compute_ss_matrix(const SparseMatrix<double> &inviscid_ss_matrix,
                                   SparseMatrix<double> &ss_matrix);

   private:
      void compute_diffusion_matrix();

      SparsityPattern      sparsity_pattern;
      SparseMatrix<double> diffusion_matrix;

      ConstraintMatrix constraints;

   protected:
      Vector<double> viscosity;

      const unsigned int n_cells;
      const unsigned int dofs_per_cell;
      const DoFHandler<dim> * const dof_handler;
};

#include "Viscosity.cc"
#endif
