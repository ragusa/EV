/**
 * \file FCT.cc
 * \brief Provides the function definitions for the FCT class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 */
template <int dim>
FCT<dim>::FCT(const RunParameters & run_parameters_,
              const DoFHandler<dim> & dof_handler_,
              const FESystem<dim> & fe_,
              const std::map<unsigned int, double> & dirichlet_values_)
  : run_parameters(&run_parameters_),
    dirichlet_values(&dirichlet_values_),
    dof_handler(&dof_handler_),
    fe(&fe_),
    n_dofs(dof_handler_.n_dofs()),
    filter_sequence_string(run_parameters_.filter_sequence_string),
    n_filters(0),
    bounds_transient_file_index(1)
{
  // create sparsity pattern for limiter and antidiffusion matrices
  DynamicSparsityPattern dsp(n_dofs);
  DoFTools::make_sparsity_pattern(*dof_handler, dsp);
  sparsity_pattern = std::make_shared<SparsityPattern>();
  sparsity_pattern->copy_from(dsp);

  // initialize matrices with sparsity pattern
  limiter_matrix.reinit(*sparsity_pattern);
  antidiffusion_matrix.reinit(*sparsity_pattern);
  limited_antidiffusion_matrix.reinit(*sparsity_pattern);

  // resize vectors
  cumulative_antidiffusion.reinit(this->n_dofs);

  // determine whether antidiffusion should be reported
  const bool report_antidiffusion = run_parameters_.verbosity_level > 1;

  // create limiter
  switch (run_parameters_.limiter_option)
  {
    case LimiterOption::ones: // no limiter; all L_{i,j} are one
      limiter = std::make_shared<OnesLimiter<dim>>(n_dofs, report_antidiffusion);
      break;
    case LimiterOption::zeroes: // full limiter; all L_{i,j} are zero
      limiter =
        std::make_shared<ZeroesLimiter<dim>>(n_dofs, report_antidiffusion);
      break;
    case LimiterOption::zalesak: // Zalesak limiter
      limiter =
        std::make_shared<ZalesakLimiter<dim>>(n_dofs, report_antidiffusion);
      break;
    default:
      AssertThrow(false, ExcNotImplemented());
      break;
  }

  // if using multi-pass limiter, then create it and switch pointer
  if (run_parameters_.use_multipass_limiting)
  {
    // switch current limiter to another pointer
    std::shared_ptr<Limiter<dim>> unit_limiter = limiter;

    // create multi-pass limiter
    limiter = std::make_shared<MultipassLimiter<dim>>(
      n_dofs,
      unit_limiter,
      sparsity_pattern,
      run_parameters_.multipass_limiting_percent_tolerance,
      report_antidiffusion);
  }
}

/**
 * \brief Computes the row sum of each row in a matrix.
 *
 * The row sum vector \f$\mathbf{a}\f$ of the matrix \f$\mathbf{A}\f$ is
 * computed as
 * \f[
 *   a_i = \sum\limits_j A_{i,j} \,.
 * \f]
 *
 * \pre This function assumes that the vector is already sized to the matrix.
 *
 * \param[in] matrix  matrix \f$\mathbf{A}\f$
 * \param[out] row_sum_vector  vector \f$\mathbf{a}\f$ of row sums of matrix
 *             \f$\mathbf{A}\f$
 */
template <int dim>
void FCT<dim>::compute_row_sum_vector(const SparseMatrix<double> & matrix,
                                      Vector<double> & row_sum_vector) const
{
  // reset row sum vector
  row_sum_vector = 0;

  // compute row sums
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    SparseMatrix<double>::const_iterator it = matrix.begin(i);
    SparseMatrix<double>::const_iterator it_end = matrix.end(i);
    for (; it != it_end; ++it)
      row_sum_vector[i] += it->value();
  }
}

/**
 * \brief Outputs bounds to files.
 *
 * \param[in] postprocessor  post-processor
 */
template <int dim>
void FCT<dim>::output_bounds(PostProcessor<dim> & postprocessor) const
{
  // use a method from the post-processor class to output the bounds
  postprocessor.output_solution(
    get_lower_solution_bound(), 0.0, *this->dof_handler, "DMPmin");
  postprocessor.output_solution(
    get_upper_solution_bound(), 0.0, *this->dof_handler, "DMPmax");
}

/**
 * \brief Outputs FCT bounds during a transient if time step counter is at
 *        a value given by the desired transient output frequency.
 *
 * \param[in] postprocessor post-processor
 * \param[in] time current time
 */
template <int dim>
void FCT<dim>::output_bounds_transient(PostProcessor<dim> & postprocessor,
                                       const double & time)
{
  // output lower bound
  postprocessor.output_dof_transient(get_lower_solution_bound(),
                                     time,
                                     *dof_handler,
                                     "lower_fct_bound",
                                     lower_bound_component_names,
                                     bounds_transient_file_index,
                                     times_and_lower_bound_filenames,
                                     false,
                                     false);

  // output upper bound
  postprocessor.output_dof_transient(get_upper_solution_bound(),
                                     time,
                                     *dof_handler,
                                     "upper_fct_bound",
                                     upper_bound_component_names,
                                     bounds_transient_file_index,
                                     times_and_upper_bound_filenames,
                                     false,
                                     false);

  // increment transient file index for bounds
  bounds_transient_file_index++;
}

/**
 * \brief Outputs the matrix of limiting coefficients.
 */
template <int dim>
void FCT<dim>::output_limiter_matrix() const
{
  // save limiting coefficients
  std::ofstream limiter_out("output/limiter.txt");
  limiter_matrix.print_formatted(limiter_out, 10, true, 0, "0", 1);
  limiter_out.close();
}

/**
 * \brief Creates a vector of filter identifier strings
 *
 * \param[in] filter_sequence_string  string of comma-delimited filter identifiers
 *
 * \return vector of filter identifier strings
 */
template <int dim>
std::vector<std::string> FCT<dim>::create_filter_string_list(
  std::string filter_sequence_string)
{
  // separate filter sequence string into vector of filter strings
  std::vector<std::string> filter_strings;
  char delimiter = ',';
  std::stringstream ss(filter_sequence_string);
  std::string filter_string;
  while (std::getline(ss, filter_string, delimiter))
  {
    filter_strings.push_back(filter_string);
  }

  // store number of filters
  n_filters = filter_strings.size();

  // return vector of filter strings
  return filter_strings;
}
