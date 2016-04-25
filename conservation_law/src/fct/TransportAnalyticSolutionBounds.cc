/**
 * \file TransportAnalyticSolutionBounds.cc
 * \brief Provides the function definitions for the
 * TransportAnalyticSolutionBounds class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] problem_parameters_  problem parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] cell_quadrature_  cell quadrature
 * \param[in] n_dofs_  number of degrees of freedom
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 */
template <int dim>
TransportAnalyticSolutionBounds<dim>::TransportAnalyticSolutionBounds(
  TransportProblemParameters<dim> & problem_parameters_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const QGauss<dim> & cell_quadrature_,
  const std::map<unsigned int, double> & dirichlet_values_)
  : DoFBounds<dim>(dof_handler_, fe_),
    cross_section_bounds(problem_parameters_.cross_section_function,
                         false,
                         dof_handler_,
                         fe_,
                         cell_quadrature_),
    source_bounds(problem_parameters_.source_function,
                  problem_parameters_.source_is_time_dependent,
                  dof_handler_,
                  fe_,
                  cell_quadrature_),
    speed(problem_parameters_.transport_speed),
    dirichlet_values(&dirichlet_values_)
{
}

/**
 * \brief Computes and applies the analytic solution bounds.
 *
 * Note that speed is assumed to be equal to one.
 *
 * \param[in] solution  solution at which to evaluate min/max
 * \param[in] dt        time step size
 * \param[in] time_old  old time
 */
template <int dim>
void TransportAnalyticSolutionBounds<dim>::update(const Vector<double> & solution,
                                                  const double & dt,
                                                  const double & time_old)
{
  // compute v*dt
  const double vdt = speed * dt;

  // update function bounds
  cross_section_bounds.update(time_old);
  source_bounds.update(time_old);

  // compute min and max of solution in neighborhood of each DoF
  this->compute_min_max_dof_vector(solution, this->lower, this->upper);

  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    // if not a Dirichlet value
    if (dirichlet_values->find(i) == dirichlet_values->end())
    {
      // branch on minimum reaction coefficient to avoid division by zero
      double min_source_term;
      if (std::abs(cross_section_bounds.upper[i]) <
          1.0e-15) // equal to zero within precision
        min_source_term = vdt * source_bounds.lower[i];
      else
        min_source_term = source_bounds.lower[i] / cross_section_bounds.upper[i] *
          (1.0 - std::exp(-vdt * cross_section_bounds.upper[i]));

      // branch on maximum reaction coefficient to avoid division by zero
      double max_source_term;
      if (std::abs(cross_section_bounds.lower[i]) <
          1.0e-15) // equal to zero within precision
        max_source_term = vdt * source_bounds.upper[i];
      else
        max_source_term = source_bounds.upper[i] / cross_section_bounds.lower[i] *
          (1.0 - std::exp(-vdt * cross_section_bounds.lower[i]));

      // compute bounds
      this->lower[i] =
        this->lower[i] * std::exp(-vdt * cross_section_bounds.upper[i]) +
        min_source_term;
      this->upper[i] =
        this->upper[i] * std::exp(-vdt * cross_section_bounds.lower[i]) +
        max_source_term;
    }
    else
    {
      // get Dirichlet value
      const double value = dirichlet_values->at(i);
      this->lower[i] = value;
      this->upper[i] = value;
    }
  }

  /*
    for (unsigned int i = 0; i < this->n_dofs; ++i)
  std::cout << this->upper[i] << std::endl;
  std::cout << std::endl;
    for (unsigned int i = 0; i < this->n_dofs; ++i)
  std::cout << this->lower[i] << std::endl;
  */
}
