/**
 * \file ExplicitEulerFCT.h
 * \brief Provides the header for the ExplicitEulerFCT class.
 */
#ifndef ExplicitEulerFCT_h
#define ExplicitEulerFCT_h

#include "include/fct/DMPExplicitEulerFCTFilter.h"
#include "include/fct/FCT.h"

using namespace dealii;

/**
 * \brief Class for explicit Euler FCT.
 */
template <int dim>
class ExplicitEulerFCT : public FCT<dim>
{
public:
  ExplicitEulerFCT(const RunParameters & run_parameters,
                   const DoFHandler<dim> & dof_handler,
                   const FESystem<dim> & fe,
                   const SparseMatrix<double> & consistent_mass_matrix,
                   const SparseMatrix<double> & lumped_mass_matrix);

  void compute_antidiffusion_vector(
    const Vector<double> & high_order_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const Vector<double> & inviscid_ss_flux,
    const Vector<double> & ss_reaction,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix,
    const Vector<double> & ss_rhs,
    Vector<double> & antidiffusion_vector);

  bool check_bounds(const Vector<double> & new_solution) const override;

  Vector<double> get_lower_solution_bound() const override;

  Vector<double> get_upper_solution_bound() const override;

protected:
  void create_filters();

  virtual std::shared_ptr<ExplicitEulerFCTFilter<dim>> create_filter(
    const std::string & filter_string);

  void compute_antidiffusion_matrix(
    const Vector<double> & high_order_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix,
    SparseMatrix<double> & antidiffusion_matrix);

  void filter_antidiffusive_fluxes(
    const Vector<double> & old_solution,
    const double & dt,
    const Vector<double> & inviscid_ss_flux,
    const Vector<double> & ss_reaction,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const Vector<double> & ss_rhs,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix);

  /** \brief vector of FCT filters */
  std::vector<std::shared_ptr<ExplicitEulerFCTFilter<dim>>> filters;

  /** \brief consistent mass matrix \f$\mathbf{M}^C\f$ */
  const SparseMatrix<double> * const consistent_mass_matrix;

  /** \brief lumped mass matrix \f$\mathbf{M}^L\f$ */
  const SparseMatrix<double> * const lumped_mass_matrix;

  /** \brief temporary vector */
  Vector<double> tmp_vector;
};

#include "src/fct/ExplicitEulerFCT.cc"

#endif
