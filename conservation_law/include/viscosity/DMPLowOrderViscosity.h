/**
 * \file DMPLowOrderViscosity.h
 * \brief Provides the header for the DMPLowOrderViscosity class.
 */
#ifndef DMPLowOrderViscosity_h
#define DMPLowOrderViscosity_h

using namespace dealii;

/**
 * \brief Class for low-order viscosity that satisfies a discrete maximum
 *        principle. This is only valid for scalar conservation laws.
 */
template <int dim>
class DMPLowOrderViscosity : public Viscosity<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  /** \brief Alias for cell iterator map to double */
  using CellMap = std::map<Cell, double>;

  DMPLowOrderViscosity(const DoFHandler<dim> & dof_handler,
                       const SparseMatrix<double> & inviscid_matrix,
                       const unsigned int & dofs_per_cell);

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n) override;

  void reinitialize() override;

private:
  void compute_viscous_bilinear_forms();

  /** \brief Sparsity pattern for viscous bilinear forms matrix */
  SparsityPattern sparsity_pattern;

  /** \brief Matrix of graph-theoretic viscous bilinear form sums */
  SparseMatrix<double> viscous_bilinear_forms;

  /** \brief Pointer to inviscid steady-state matrix */
  const SparseMatrix<double> * const inviscid_matrix;

  /** \brief Number of degrees of freedom per cell */
  const unsigned int dofs_per_cell;
};

#include "src/viscosity/DMPLowOrderViscosity.cc"

#endif
