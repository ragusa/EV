/**
 * \file FCT.cc
 * \brief Provides the function definitions for the FCT class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] dof_handler_  degree of freedom handler
 */
template <int dim>
FCT<dim>::FCT(const RunParameters & run_parameters_,
              const DoFHandler<dim> & dof_handler_)
  : dof_handler(&dof_handler_),
    n_dofs(dof_handler_.n_dofs()),
    filter_sequence_string(run_parameters_.filter_sequence_string)
{
  // create limiter
  switch (run_parameters_.limiter_option)
  {
    case LimiterOption::ones: // no limiter; all L_{i,j} are one
      limiter = std::make_shared<OnesLimiter>(n_dofs);
      break;
    case LimiterOption::zeroes: // full limiter; all L_{i,j} are zero
      limiter = std::make_shared<ZeroesLimiter>(n_dofs);
      break;
    case LimiterOption::zalesak: // Zalesak limiter
      limiter = std::make_shared<ZalesakLimiter>(n_dofs);
      break;
    default:
      throw ExcNotImplemented();
      break;
  }

  // create sparsity pattern for limiter and antidiffusion matrices
  DynamicSparsityPattern dsp(n_dofs);
  DoFTools::make_sparsity_pattern(*dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  // initialize matrices with sparsity pattern
  limiter_matrix.reinit(sparsity_pattern);
  antidiffusion_matrix.reinit(sparsity_pattern);
  limited_antidiffusion_matrix.reinit(sparsity_pattern);

  // resize vectors
  cumulative_antidiffusion.reinit(this->n_dofs);
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
