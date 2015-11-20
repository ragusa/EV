/**
 * \file Entropy.cc
 * \brief Provides the function definitions for the Entropy class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] use_max_entropy_deviation_normalization_ flag to use max
 *            entropy deviation from average as entropy normalization
 * \param[in] domain_volume_ volume of domain
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] fe_ finite element system
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] face_quadrature_ face quadrature
 */
template <int dim>
Entropy<dim>::Entropy(const bool & use_max_entropy_deviation_normalization_,
                      const double & domain_volume_,
                      const DoFHandler<dim> & dof_handler_,
                      const FESystem<dim> & fe_,
                      const QGauss<dim> & cell_quadrature_,
                      const QGauss<dim - 1> & face_quadrature_)
  : dof_handler(&dof_handler_),
    fe(&fe_),
    cell_quadrature(&cell_quadrature_),
    face_quadrature(&face_quadrature_),
    n_q_points_cell(cell_quadrature_.size()),
    n_q_points_face(face_quadrature_.size()),
    use_max_entropy_deviation_normalization(
      use_max_entropy_deviation_normalization_),
    domain_volume(domain_volume_)
{
}

/**
 * \brief Reinitializes group FE values (if any) and computes
 *        solution-dependent quantities for normalization of entropy
 *        viscosity.
 *
 * \param[in] solution solution vector
 */
template <int dim>
void Entropy<dim>::reinitialize(const Vector<double> & solution)
{
  // reinitialize group FE values (if any)
  reinitialize_group_fe_values(solution);

  // compute average entropy and max entropy deviation if needed
  if (use_max_entropy_deviation_normalization)
  {
    // compute max entropy deviation
    compute_average_entropy(solution);
    compute_max_entropy_deviation(solution);
  }
}

/**
 * \brief Reinitializes group FE values, if any.
 *
 * \param[in] solution solution vector
 */
template <int dim>
void Entropy<dim>::reinitialize_group_fe_values(const Vector<double> &)
{
}

/**
 * Computes the average of the entropy over the domain.
 *
 * \param[in] solution solution with which to calculate entropy for the
 *            normalization constant
 */
template <int dim>
void Entropy<dim>::compute_average_entropy(const Vector<double> & solution)
{
  FEValues<dim> fe_values(
    *fe, *cell_quadrature, update_values | update_JxW_values);

  double domain_integral_entropy = 0.0;

  // loop over cells
  Cell cell = dof_handler->begin_active(), endc = dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // compute entropy
    auto entropy = compute_entropy(solution, fe_values);

    // add contribution of quadrature point to entropy integral
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      domain_integral_entropy += entropy[q] * fe_values.JxW(q);
  }

  // finish computing average entropy
  entropy_average = domain_integral_entropy / domain_volume;
}

/**
 * \brief Computes the max entropy deviation from the average of each quadrature
 *        point in the domain.
 *
 * \param[in] solution solution with which to calculate entropy for the
 *            normalization constant
 */
template <int dim>
void Entropy<dim>::compute_max_entropy_deviation(const Vector<double> & solution)
{
  FEValues<dim> fe_values(*fe, *cell_quadrature, update_values);

  // initialize max entropy deviation to zero
  max_entropy_deviation = 0.0;

  // loop over cells
  Cell cell = dof_handler->begin_active(), endc = dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // compute entropy
    auto entropy = compute_entropy(solution, fe_values);

    // loop over quadrature points
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      max_entropy_deviation =
        std::max(max_entropy_deviation, std::abs(entropy[q] - entropy_average));
  }
}

/**
 * \brief Returns the maximum entropy deviation from the average as a uniform
 *        vector, to be used as the entropy normalization constant.
 *
 * \param[in] solution solution with which to calculate entropy for the
 *            normalization constant
 * \param[in] cell cell iterator
 *
 * \return vector of maximum entropy deviation
 */
template <int dim>
std::vector<double> Entropy<dim>::compute_max_entropy_deviation_normalization(
  const Vector<double> &, const Cell &) const
{
  std::vector<double> normalization(n_q_points_cell, max_entropy_deviation);

  return normalization;
}
