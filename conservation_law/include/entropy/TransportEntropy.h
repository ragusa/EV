/**
 * \file TransportEntropy.h
 * \brief Provides the header for the TransportEntropy class.
 */
#ifndef TransportEntropy_h
#define TransportEntropy_h

#include "include/parameters/TransportProblemParameters.h"
#include "include/entropy/Entropy.h"

using namespace dealii;

/**
 * \brief Class for entropy for the scalar transport equation.
 *
 * For now, the following entropy function is chosen:
 * \f[
 *   \eta(u) = \frac{1}{2} u^2 \,.
 * \f]
 * The entropy residual is the following:
 * \f[
 *   \mathcal{R}_K^n =
 *     \frac{\eta(\tilde{u}^n_q) - \eta(\tilde{u}^{n-1}_q)}{\Delta t^{n-1}}
 *     + \eta'(\tilde{u}_q)v(\mathbf{\Omega}\cdot\nabla \tilde{u}^n_q
 *       + \sigma_q\tilde{u}^n_q - q^n_q) \,.
 * \f]
 */
template <int dim>
class TransportEntropy : public Entropy<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Entropy<dim>::Cell;

  TransportEntropy(const double & domain_volume,
                   const DoFHandler<dim> & dof_handler,
                   const FESystem<dim> & fe,
                   const QGauss<dim> & cell_quadrature,
                   const QGauss<dim - 1> & face_quadrature,
                   const TransportProblemParameters<dim> & problem_parameters);

  std::vector<double> compute_entropy_residual(
    const Vector<double> & new_solution,
    const Vector<double> & old_solution,
    const double & dt,
    const Cell & cell) override;

  double compute_max_entropy_jump(const Vector<double> & solution,
                                  const Cell & cell) override;

private:
  std::vector<double> compute_entropy(
    const Vector<double> & solution,
    const FEValuesBase<dim> & fe_values) const override;

  std::vector<double> compute_entropy_scalar(
    const std::vector<double> & solution_local) const noexcept;

  std::vector<double> compute_entropy_derivative(
    const std::vector<double> & solution_local) const noexcept;

  /** \brief Transport speed \f$v\f$ */
  const double transport_speed;

  /** \brief Transport direction \f$\mathbf{\Omega}\f$ */
  const Tensor<1, dim> transport_direction;

  /** \brief Function parser of cross section \f$\sigma\f$ */
  const FunctionParser<dim> * const cross_section_function;

  /** \brief Function parser of source \f$q\f$ */
  const FunctionParser<dim> * const source_function;
};

#include "src/entropy/TransportEntropy.cc"

#endif
