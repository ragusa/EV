/**
 * \file TransportUpwindAnalyticSolutionBounds.cc
 * \brief Provides the function definitions for the
 * TransportUpwindAnalyticSolutionBounds class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] problem_parameters_  problem parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] n_samples_  number of sampling points for min/max determination
 */
template <int dim>
TransportUpwindAnalyticSolutionBounds<dim>::TransportUpwindAnalyticSolutionBounds(
  TransportProblemParameters<dim> & problem_parameters_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const unsigned int & n_samples_)
  : DoFBounds<dim>(dof_handler_, fe_),
    cross_section_function(&problem_parameters_.cross_section_function),
    source_function(&problem_parameters_.source_function),
    incoming(problem_parameters_.incoming),
    direction(problem_parameters_.transport_direction),
    boundary_distance(problem_parameters_.boundary_distance),
    n_samples(n_samples_)
{
  // get vector of positions of each DoF
  x.resize(this->n_dofs);
  DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), *this->dof_handler, x);

  // resize sampling vectors and compute sample fractions
  x_sample.resize(n_samples);
  alpha_sample.resize(n_samples);
  for (unsigned int j = 0; j < n_samples; ++j)
    alpha_sample[j] = 1.0 / (n_samples + 1) * j;
  sigma_samples.resize(n_samples);
  source_samples.resize(n_samples);
}

/**
 * \brief Computes and applies the analytic solution bounds.
 *
 * These solution bounds are computed as
 * \f[
 *   W_i^\pm = u_{r,i}^\pm + u_{q,i}^\pm ,
 * \f]
 * \f[
 *   u_{r,i}^- = u_{ref,i} e^{-s_{ref,i}\sigma_i^+} , \quad
 *   u_{r,i}^+ = u_{ref,i} e^{-s_{ref,i}\sigma_i^-} ,
 * \f]
 * \f[
 *   u_{q,i}^- = \left\{\begin{array}{l l}
 *     \frac{q_i^-}{\sigma_i^+}\left(1 - e^{-s_{ref,i}\sigma_i^+}\right)
 *       & \sigma_i^+ \ne 0\\
 *     s_{ref,i}q_i^-
 *       & \sigma_i^+ = 0\\
 *   \end{array}\right.
 * \f]
 * \f[
 *   u_{q,i}^+ = \left\{\begin{array}{l l}
 *     \frac{q_i^+}{\sigma_i^-}\left(1 - e^{-s_{ref,i}\sigma_i^-}\right)
 *       & \sigma_i^- \ne 0\\
 *     s_{ref,i}q_i^+
 *       & \sigma_i^- = 0\\
 *   \end{array}\right.
 * \f]
 * \f[
 *   \sigma_i^- = \min\limits_j \sigma(\mathbf{x}_j^{sample}) , \quad
 *   \sigma_i^+ = \max\limits_j \sigma(\mathbf{x}_j^{sample}) ,
 * \f]
 * \f[
 *   q_i^- = \min\limits_j q(\mathbf{x}_j^{sample}) , \quad
 *   q_i^+ = \max\limits_j q(\mathbf{x}_j^{sample}) ,
 * \f]
 * \f[
 *   \mathbf{x}_j^{sample} = \mathbf{x}_i - \alpha_j s_{ref,i} \mathbf{\Omega} ,
 * \quad
 *   \alpha_j = \frac{1}{N_{sample}+1}j ,
 * \f]
 * \f[
 *   s_{ref,i} = \min(s, s_{b,i}) ,
 * \f]
 * \f[
 *   u_{ref,i} = u(\mathbf{x}_i - s_{ref,i}\mathbf{\Omega}) ,
 * \f]
 * where \f$s_{b,i}\f$ is the distance from \f$\mathbf{x}_i\f$ to the incoming
 * boundary along \f$\mathbf{\Omega}\f$.
 *
 * \param[in] solution  solution at which to evaluate min/max
 * \param[in] s         path length (\f$v\Delta t\f$ for transient,
 *            \f$\Delta x\f$ for steady-state
 * \param[in] time_old  old time
 */
template <int dim>
void TransportUpwindAnalyticSolutionBounds<dim>::update(
  const Vector<double> & solution, const double & s, const double & time_old)
{
  // set time for source
  source_function->set_time(time_old);

  // loop over DoFs
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    // compute distance to boundary
    const double sb = boundary_distance->compute(x[i]);

    // determine distance to reference point
    const double sref = std::min(s, sb);

    // compute reference position
    Point<dim> xref = x[i] - sref * direction;

    // determine reference solution
    double uref;
    if (s >= sb) // if the distance reaches the boundary, use incoming BC value
      uref = incoming;
    else
      // compute the reference value from solution
      uref = VectorTools::point_value(*this->dof_handler, solution, xref);

    // create sample points
    for (unsigned int j = 0; j < n_samples; ++j)
      x_sample[j] = x[i] - alpha_sample[j] * sref * direction;

    // sample cross section and source
    cross_section_function->value_list(x_sample, sigma_samples);
    source_function->value_list(x_sample, source_samples);

    // determine min and max cross section and source
    const double sigma_min =
      *std::min_element(sigma_samples.begin(), sigma_samples.end());
    const double sigma_max =
      *std::max_element(sigma_samples.begin(), sigma_samples.end());
    const double source_min =
      *std::min_element(source_samples.begin(), source_samples.end());
    const double source_max =
      *std::max_element(source_samples.begin(), source_samples.end());

    // compute bounds on component of solution associated with reference value
    const double ur_min = uref * std::exp(-sref * sigma_max);
    const double ur_max = uref * std::exp(-sref * sigma_min);

    // compute bounds on component of solution associated with source
    double uq_min;
    if (sigma_max > 1.0e-15)
      uq_min = source_min / sigma_max * (1.0 - std::exp(-sref * sigma_max));
    else
      uq_min = sref * source_min;

    double uq_max;
    if (sigma_min > 1.0e-15)
      uq_max = source_max / sigma_min * (1.0 - std::exp(-sref * sigma_min));
    else
      uq_max = sref * source_max;

    // put components together
    this->lower[i] = ur_min + uq_min;
    this->upper[i] = ur_max + uq_max;
  }
}
