/**
 * \file MultipassLimiter.h
 * \brief Provides the header for the MultipassLimiter class.
 */
#ifndef MultipassLimiter_h
#define MultipassLimiter_h

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "include/fct/Limiter.h"

using namespace dealii;

/**
 * \brief Class for multi-pass limiter.
 *
 * This limiter allows multiple sweeps through any other limiter to allow
 * more antidiffusion to be kept.
 */
template <int dim>
class MultipassLimiter : public Limiter<dim>
{
public:
  MultipassLimiter(const unsigned int & n_dofs,
                   std::shared_ptr<Limiter<dim>> limiter,
                   const std::shared_ptr<SparsityPattern> sparsity_pattern,
                   const double & percent_tolerance = 0.01,
                   const bool & report_antidiffusion = false);

  void compute_limiter_matrix(
    const SparseMatrix<double> & antidiffusion_matrix,
    const DoFBounds<dim> & antidiffusion_bounds,
    const Vector<double> & cumulative_antidiffusion_vector,
    SparseMatrix<double> & limiter_matrix) override;

  void apply_limiter_matrix(const SparseMatrix<double> & limiter_matrix,
                            SparseMatrix<double> & antidiffusion_matrix) const;

protected:
  void update_remainder_matrix_and_accepted_antidiffusion(
    const SparseMatrix<double> & limiter_matrix);

  void compute_total_limiter_matrix(
    const SparseMatrix<double> & antidiffusion_matrix,
    const SparseMatrix<double> & remainder_antidiffusion_matrix,
    SparseMatrix<double> & limiter_matrix) const;

  /** \brief limiter for which multiple passes are to be made */
  std::shared_ptr<Limiter<dim>> limiter;

  /** \brief sparsity pattern */
  const std::shared_ptr<SparsityPattern> sparsity_pattern;

  /** \brief percent tolerance for ending passes */
  const double percent_tolerance;

  /** \brief remainder antidiffusion matrix \f$\Delta\mathbf{P}\f$ */
  SparseMatrix<double> remainder_matrix;

  /** \brief accepted antidiffusion \f$\bar{\mathbf{p}}\f$ */
  Vector<double> accepted_antidiffusion;
};

#include "src/fct/MultipassLimiter.cc"

#endif
