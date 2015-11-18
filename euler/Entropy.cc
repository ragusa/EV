/**
 * \file Entropy.cc
 * \brief Provides the function definitions for the Entropy class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
Entropy<dim>::Entropy()
  : need_to_compute_average_entropy()
{
}

/**
 * \brief Reinitializes group FE values (if any) and computes
 *        solution-dependent quantities for normalization of entropy
 *        viscosity.
 *
 * \param[in] solution solution vector
 */
void Entropy<dim>::reinitialize(const Vector<double> & solution)
{
  // reinitialize group FE values (if any)
  reinitialize_group_fe_values(solution);

  // compute average entropy if needed for normalization of entropy viscosity
  if (need_to_compute_average_entropy)
    entropy_average = compute_average_entropy(solution);
}

/**
 * Computes the average of the entropy over the domain.
 *
 * \param[in] solution solution with which to calculate entropy for the
 *            normalization constant
 *
 * \return average entropy in domain
 */
template <int dim>
double Entropy<dim>::compute_average_entropy(
  const Vector<double> & solution) const
{
  FEValues<dim> fe_values(
    *fe, cell_quadrature, update_values | update_JxW_values);

  Vector<double> entropy(n_q_points_cell);

  double domain_integral_entropy = 0.0;

  // loop over cells
  Cell cell = dof_handler.begin_active(), endc = dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // compute entropy
    compute_entropy(solution, fe_values, entropy);

    // loop over quadrature points
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      // add contribution of quadrature point to entropy integral
      domain_integral_entropy += entropy[q] * fe_values.JxW(q);
    }
  }
  // domain-averaged entropy_values
  const double domain_averaged_entropy = domain_integral_entropy / domain_volume;

  return domain_averaged_entropy;
}

/**
 * \brief Computes the entropy deviation from the average at each quadrature
 *        point in a cell.
 *
 * \param[in] solution solution with which to calculate entropy for the
 *            normalization constant
 * \param[in] entropy_average average entropy in domain
 *
 * \return vector of entropy deviations at each quadrature point in cell
 */
template <int dim>
std::vector<double> Entropy<dim>::compute_entropy_normalization(
  const Vector<double> & solution,
  const double & entropy_average,
  const Cell & cell) const
{
  FEValues<dim> fe_values(fe, cell_quadrature, update_values);

  Vector<double> entropy(n_q_points_cell);

  std::vector<double> normalization_constant(n_q_points_cell);

  // reinitialize FE values
  fe_values.reinit(cell);

  // compute entropy
  compute_entropy(solution, fe_values, entropy);

  // loop over quadrature points
  for (unsigned int q = 0; q < n_q_points_cell; ++q)
    normalization_constant[q] = std::abs(entropy[q] - entropy_average);

  return normalization_constant;
}
