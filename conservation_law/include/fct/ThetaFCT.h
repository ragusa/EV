/**
 * \file ThetaFCT.h
 * \brief Provides the header for the ThetaFCT class.
 */
#ifndef ThetaFCT_h
#define ThetaFCT_h

#include "include/fct/DMPThetaFCTFilter.h"
#include "include/fct/FCT.h"

using namespace dealii;

/**
 * \brief Class for theta FCT.
 */
template <int dim>
class ThetaFCT : public FCT<dim>
{
public:
  ThetaFCT(const RunParameters & run_parameters,
           const DoFHandler<dim> & dof_handler,
           const FESystem<dim> & fe,
           const SparseMatrix<double> & consistent_mass_matrix,
           const SparseMatrix<double> & lumped_mass_matrix);

  void compute_antidiffusion_matrix(
    const Vector<double> & high_order_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_diffusion_matrix_old,
    const SparseMatrix<double> & low_order_diffusion_matrix_new,
    const SparseMatrix<double> & high_order_diffusion_matrix_old,
    const SparseMatrix<double> & high_order_diffusion_matrix_new);

  void compute_antidiffusion_vector(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs_new,
    const Vector<double> & ss_rhs_old,
    Vector<double> & antidiffusion_vector);

  bool check_bounds(const Vector<double> & new_solution) const override;

  Vector<double> get_lower_solution_bound() const override;

  Vector<double> get_upper_solution_bound() const override;

protected:
  void create_filters();

  virtual std::shared_ptr<ThetaFCTFilter<dim>> create_filter(
    const std::string & filter_string);

  void filter_antidiffusive_fluxes(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const SparseMatrix<double> & low_order_ss_matrix,
    const Vector<double> & ss_rhs_new,
    const Vector<double> & ss_rhs_old,
    const Vector<double> & cumulative_antidiffusion,
    SparseMatrix<double> & limiter_matrix,
    SparseMatrix<double> & antidiffusion_matrix);

  /** \brief vector of FCT filters */
  std::vector<std::shared_ptr<ThetaFCTFilter<dim>>> filters;

  /** \brief consistent mass matrix \f$\mathbf{M}^C\f$ */
  const SparseMatrix<double> * const consistent_mass_matrix;

  /** \brief lumped mass matrix \f$\mathbf{M}^L\f$ */
  const SparseMatrix<double> * const lumped_mass_matrix;

  /** \brief theta parameter \f$\theta\f$ */
  const double theta;

  /** \brief flag to use cumulative antidiffusion algorithm */
  const bool use_cumulative_antidiffusion_algorithm;

  /** \brief temporary vector */
  Vector<double> tmp_vector;
};

#include "src/fct/ThetaFCT.cc"

#endif
