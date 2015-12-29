/**
 * \file FCT.cc
 * \brief Provides the function definitions for the FCT class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
FCT<dim>::FCT(const ConservationLawParameters<dim> & parameters_,
              const DoFHandler<dim> & dof_handler_,
              const Triangulation<dim> & triangulation_,
              const SparseMatrix<double> & lumped_mass_matrix_,
              const SparseMatrix<double> & consistent_mass_matrix_,
              const std::shared_ptr<StarState<dim>> & star_state_,
              const LinearSolver<dim> & linear_solver_,
              const SparsityPattern & sparsity_pattern_,
              const std::vector<unsigned int> & dirichlet_nodes_,
              const unsigned int & n_components_,
              const unsigned int & dofs_per_cell_)
  : fe_scalar(1),
    dof_handler(&dof_handler_),
    dof_handler_scalar(triangulation_),
    lumped_mass_matrix(&lumped_mass_matrix_),
    consistent_mass_matrix(&consistent_mass_matrix_),
    star_state(star_state_),
    linear_solver(linear_solver_),
    sparsity(&sparsity_pattern_),
    dirichlet_nodes(dirichlet_nodes_),
    n_dofs(dof_handler_.n_dofs()),
    n_components(n_components_),
    n_dofs_scalar(n_dofs / n_components),
    dofs_per_cell(dofs_per_cell_),
    dofs_per_cell_per_component(dofs_per_cell_ / n_components_),
    antidiffusion_type(parameters_.antidiffusion_type),
    synchronization_type(parameters_.fct_synchronization_type),
    fct_variables_type(parameters_.fct_variables_type),
    use_star_states_in_fct_bounds(parameters_.use_star_states_in_fct_bounds),
    DMP_satisfied_at_all_steps(true)
{
  // distribute degrees of freedom in scalar degree of freedom handler
  dof_handler_scalar.clear();
  dof_handler_scalar.distribute_dofs(fe_scalar);

  // create scalar sparsity pattern
  DynamicSparsityPattern dsp(n_dofs_scalar);
  DoFTools::make_sparsity_pattern(dof_handler_scalar, dsp);
  scalar_sparsity_pattern.copy_from(dsp);

  // allocate memory for vectors
  tmp_vector.reinit(n_dofs);
  system_rhs.reinit(n_dofs);
  solution_min.reinit(n_dofs);
  solution_max.reinit(n_dofs);
  flux_correction_vector.reinit(n_dofs);
  Q_minus.reinit(n_dofs);
  Q_plus.reinit(n_dofs);
  R_minus.reinit(n_dofs);
  R_plus.reinit(n_dofs);

  // initialize sparse matrices
  system_matrix.reinit(sparsity_pattern_);
  limiter_matrix.reinit(sparsity_pattern_);
  flux_correction_matrix.reinit(sparsity_pattern_);
}

/**
 * \brief Solves the FCT system.
 *
 * The FCT system is
 * \f[
 *   \mathrm{M}^L\mathrm{U}^{n+1} = \mathrm{M}^L\mathrm{U}^n + \Delta t\left(
 *     \mathrm{b}^n - \mathbf{C}\mathbf{F}(\mathrm{U}^n)
 *     - \mathrm{D}^L(\mathrm{U}^n)\mathrm{U}^n - \bar{\mathrm{p}}\right) \,.
 * \f]
 *
 * \param[in,out] new_solution high-order solution when passed in, FCT solution
 *                when passed out
 * \param[in] old_solution old solution
 */
template <int dim>
void FCT<dim>::solve_fct_system(
  Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const Vector<double> & ss_flux,
  const Vector<double> & ss_reaction,
  const Vector<double> & ss_rhs,
  const double & dt,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const SparseMatrix<double> & high_order_diffusion_matrix)
{
  // compute flux corrections
  compute_flux_corrections(new_solution,
                           old_solution,
                           dt,
                           low_order_diffusion_matrix,
                           high_order_diffusion_matrix);

  // compute max principle min and max values
  compute_bounds(old_solution, ss_reaction, ss_rhs, dt);

  // compute limited flux correction sum and add it to rhs
  switch (antidiffusion_type)
  {
    case AntidiffusionType::limited:
      compute_limiting_coefficients_zalesak(
        old_solution, ss_flux, ss_rhs, low_order_diffusion_matrix, dt);
      switch (synchronization_type)
      {
        case FCTSynchronizationType::none:
          break;
        case FCTSynchronizationType::min:
          synchronize_min();
          break;
        case FCTSynchronizationType::compound:
          Assert(false, ExcNotImplemented());
          break;
        default:
          Assert(false, ExcNotImplemented());
          break;
      }
      compute_limited_flux_correction_vector();
      break;
    case AntidiffusionType::full:
      compute_full_flux_correction();
      break;
    case AntidiffusionType::none:
      flux_correction_vector = 0;
      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // check limited flux correction sum if requested
  if (true)
    check_limited_flux_correction_sum(flux_correction_vector);

  // form rhs: system_rhs = M*u_old + dt*(ss_rhs - ss_flux - D*u_old + f)
  system_rhs = 0;
  lumped_mass_matrix->vmult(tmp_vector, old_solution);
  system_rhs.add(1.0, tmp_vector);
  system_rhs.add(dt, ss_rhs);
  system_rhs.add(-dt, ss_flux);
  low_order_diffusion_matrix.vmult(tmp_vector, old_solution);
  system_rhs.add(-dt, tmp_vector);
  system_rhs.add(dt, flux_correction_vector);

  // solve the linear system M*u_new = system_rhs
  system_matrix.copy_from(*lumped_mass_matrix);
  linear_solver.solve(system_matrix, new_solution, system_rhs);

  /*
  // check that local discrete maximum principle is satisfied at all time steps
  bool DMP_satisfied_this_step =
    check_max_principle(new_solution, low_order_ss_matrix, dt);
  DMP_satisfied_at_all_steps =
    DMP_satisfied_at_all_steps and DMP_satisfied_this_step;
  */
}

/**
 * \brief Computes bounds to be imposed on the FCT solution, \f$W_i^-\f$ and
 *        \f$W_i^+\f$.
 *
 * The imposed bounds are computed as
 * \f[
 *   W_i^-(\mathrm{U}^n) = U_{min,i}^n\left(1
 *     - \frac{\Delta t}{M^L_{i,i}}\sigma_i\right)
 *     + \frac{\Delta t}{M^L_{i,i}}b_i^n \,,
 * \f]
 * \f[
 *   W_i^+(\mathrm{U}^n) = U_{max,i}^n\left(1
 *     - \frac{\Delta t}{M^L_{i,i}}\sigma_i\right)
 *     + \frac{\Delta t}{M^L_{i,i}}b_i^n \,,
 * \f]
 * where \f$U_{min,i}^n\equiv\min\limits_j U_j^n\f$ and
 * \f$U_{max,i}^n\equiv\max\limits_j U_j^n\f$.
 */
template <int dim>
void FCT<dim>::compute_bounds(const Vector<double> & old_solution,
                              const Vector<double> & ss_reaction,
                              const Vector<double> & ss_rhs,
                              const double & dt)
{
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    solution_max(i) = old_solution(i);
    solution_min(i) = old_solution(i);
  }

  // loop over cells
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler->begin_active(),
                                                 endc = dof_handler->end();
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // loop over solution components
    for (unsigned int m = 0; m < n_components; ++m)
    {
      // find min and max values on cell for this component
      double max_cell = old_solution(local_dof_indices[m]);
      double min_cell = old_solution(local_dof_indices[m]);
      for (unsigned int j = 0; j < dofs_per_cell_per_component; ++j)
      {
        // consider old states
        const double value_j =
          old_solution(local_dof_indices[j * n_components + m]);
        max_cell = std::max(max_cell, value_j);
        min_cell = std::min(min_cell, value_j);
      }

      // update the max and min values of neighborhood of each dof
      for (unsigned int j = 0; j < dofs_per_cell_per_component; ++j)
      {
        unsigned int i = local_dof_indices[j * n_components + m]; // global index
        solution_max(i) = std::max(solution_max(i), max_cell);
        solution_min(i) = std::min(solution_min(i), min_cell);
      }
    }
  }

  // consider star states in bounds if specified
  if (use_star_states_in_fct_bounds)
  {
    for (unsigned int i = 0; i < n_dofs; ++i)
    {
      // iterate over sparsity pattern to get degree of freedom indices
      // in the support of i
      SparsityPattern::iterator it = sparsity->begin(i);
      SparsityPattern::iterator it_end = sparsity->end(i);
      for (; it != it_end; ++it)
      {
        // get column index
        const unsigned int j = it->column();

        if (j != i)
        {
          // determine if j corresponds to the same component as i
          if (i % n_components == j % n_components)
          {
            // get star state value associated with i and j
            const double star_state_value = star_state->get_star_state(i, j);
            solution_max(i) = std::max(solution_max(i), star_state_value);
            solution_min(i) = std::min(solution_min(i), star_state_value);
          }
        }
      }
    }
  }

  // At this point, the min/max values of the old solution in the support
  // of test function i are stored in solution_min(i) and solution_max(i).
  // Now these values are multiplied by (1-dt/m(i))*sigma(i) and
  // added to dt/m(i)*b(i).

  // compute the upper and lower bounds for the maximum principle
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    solution_max(i) = solution_max(i) *
        (1.0 - dt / (*lumped_mass_matrix)(i, i) * ss_reaction(i)) +
      dt / (*lumped_mass_matrix)(i, i) * ss_rhs(i);

    solution_min(i) = solution_min(i) *
        (1.0 - dt / (*lumped_mass_matrix)(i, i) * ss_reaction(i)) +
      dt / (*lumped_mass_matrix)(i, i) * ss_rhs(i);
  }
}

/**
 * \brief Assembles the flux correction matrix \f$\mathrm{P}\f$.
 *
 * \param [in] dt current time step size
 */
template <int dim>
void FCT<dim>::compute_flux_corrections(
  const Vector<double> & high_order_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const SparseMatrix<double> & high_order_diffusion_matrix)
{
  // reset flux correction matrix to zero
  flux_correction_matrix = 0;

  // compute time derivative of high-order solution
  Vector<double> & dUdt = tmp_vector;
  dUdt = 0;
  dUdt.add(1.0 / dt, high_order_solution, -1.0 / dt, old_solution);

  // iterate over sparse matrix entries
  SparseMatrix<double>::const_iterator it_mass = consistent_mass_matrix->begin();
  SparseMatrix<double>::const_iterator it_low =
    low_order_diffusion_matrix.begin();
  SparseMatrix<double>::const_iterator it_high =
    high_order_diffusion_matrix.begin();
  SparseMatrix<double>::iterator it_flux = flux_correction_matrix.begin();
  SparseMatrix<double>::iterator it_end = flux_correction_matrix.end();

  for (; it_flux != it_end; ++it_flux, ++it_mass, ++it_low, ++it_high)
  {
    // get row and column indices
    const unsigned int i = it_flux->row();
    const unsigned int j = it_flux->column();

    // get values
    const double Mij = it_mass->value();
    const double DLij = it_low->value();
    const double DHij = it_high->value();

    // compute flux correction entry
    const double Fij = -Mij * (dUdt(j) - dUdt(i)) +
      (DLij - DHij) * (old_solution(j) - old_solution(i));

    // store value
    it_flux->value() = Fij;
  }
}

/**
 * \brief Computes the limiting coefficient vectors \f$R^+\f$ and \f$R^-\f$,
 *        used for computing the high-order solution from the low-order
 *        solution.
 */
template <int dim>
void FCT<dim>::compute_limiting_coefficients_zalesak(
  const Vector<double> & old_solution,
  const Vector<double> & ss_flux,
  const Vector<double> & ss_rhs,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const double & dt)
{
  // reset limiter matrix
  limiter_matrix = 0;

  // start computing Q+
  Q_plus = 0;
  lumped_mass_matrix->vmult(tmp_vector, old_solution);
  Q_plus.add(-1.0 / dt, tmp_vector);
  Q_plus.add(1.0, ss_flux);
  low_order_diffusion_matrix.vmult(tmp_vector, old_solution);
  Q_plus.add(1.0, tmp_vector);
  Q_plus.add(-1.0, ss_rhs);

  // copy current contents of Q+ as these components are identical
  Q_minus = Q_plus;

  // finish computing Q+ and Q-
  lumped_mass_matrix->vmult(tmp_vector, solution_max);
  Q_plus.add(1.0 / dt, tmp_vector);
  lumped_mass_matrix->vmult(tmp_vector, solution_min);
  Q_minus.add(1.0 / dt, tmp_vector);

  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // for Dirichlet nodes, set R_minus and R_plus to 1 because no limiting is
    // needed
    if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) !=
        dirichlet_nodes.end())
    {
      R_minus(i) = 1.0;
      R_plus(i) = 1.0;
    }
    else
    {
      // get nonzero entries in row i of flux correction matrix
      std::vector<double> row_values;
      std::vector<unsigned int> row_indices;
      unsigned int n_col;
      get_matrix_row(flux_correction_matrix, i, row_values, row_indices, n_col);

      double P_plus_i = 0.0;
      double P_minus_i = 0.0;
      // compute P_plus, P_minus
      for (unsigned int k = 0; k < n_col; ++k)
      {
        // get value of nonzero entry k
        double Fij = row_values[k];

        P_plus_i += std::max(0.0, Fij);
        P_minus_i += std::min(0.0, Fij);
      }

      // compute R_plus(i)
      if (P_plus_i != 0.0)
        R_plus(i) = std::min(1.0, Q_plus(i) / P_plus_i);
      else
        R_plus(i) = 1.0;

      // compute R_minus(i)
      if (P_minus_i != 0.0)
        R_minus(i) = std::min(1.0, Q_minus(i) / P_minus_i);
      else
        R_minus(i) = 1.0;
    }
  }

  // compute limited flux correction sum
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // get values and indices of nonzero entries in row i of coefficient matrix A
    std::vector<double> row_values;
    std::vector<unsigned int> row_indices;
    unsigned int n_col;
    get_matrix_row(flux_correction_matrix, i, row_values, row_indices, n_col);

    // perform flux correction for dof i
    // Note that flux correction sum is sum_{j in I(S_i)} of Lij*Fij.
    for (unsigned int k = 0; k < n_col; ++k)
    {
      unsigned int j = row_indices[k];
      double Fij = row_values[k];

      // compute limiting coefficient Lij
      double Lij;
      if (Fij >= 0.0)
        Lij = std::min(R_plus(i), R_minus(j));
      else
        Lij = std::min(R_minus(i), R_plus(j));
      limiter_matrix.set(i, j, Lij);
    }
  }
}

/**
 * \brief Computes limited flux correction vector.
 */
template <int dim>
void FCT<dim>::compute_limited_flux_correction_vector()
{
  // reset vector
  flux_correction_vector = 0;

  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // iterate over nonzero columns of matrix
    SparseMatrix<double>::const_iterator it_lim = limiter_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_lim_end = limiter_matrix.end(i);
    SparseMatrix<double>::const_iterator it_flux =
      flux_correction_matrix.begin(i);
    for (; it_lim != it_lim_end; ++it_lim, ++it_flux)
    {
      const double limiting_coefficient = it_lim->value();
      const double flux = it_flux->value();
      flux_correction_vector(i) += limiting_coefficient * flux;
    }
  }
}

/**
 * \brief Computes the full flux correction vector (without any limitation)
 */
template <int dim>
void FCT<dim>::compute_full_flux_correction()
{
  // reset vector
  flux_correction_vector = 0;

  // iterate over sparse matrix entries
  SparseMatrix<double>::iterator it = flux_correction_matrix.begin();
  SparseMatrix<double>::iterator it_end = flux_correction_matrix.end();

  for (; it != it_end; ++it)
  {
    // get row index
    const unsigned int i = it->row();

    // add entry to flux correction vector
    flux_correction_vector(i) += it->value();
  }
}

/**
 * \brief Gets the values and indices of nonzero elements in a row of a sparse
 *        matrix.
 *  \param [in] matrix sparse matrix whose row will be retrieved
 *  \param [in] i index of row to be retrieved
 *  \param [out] row_values vector of values of nonzero entries of row i
 *  \param [out] row_indices vector of indices of nonzero entries of row i
 *  \param [out] n_col number of nonzero entries of row i
 */
template <int dim>
void FCT<dim>::get_matrix_row(const SparseMatrix<double> & matrix,
                              const unsigned int & i,
                              std::vector<double> & row_values,
                              std::vector<unsigned int> & row_indices,
                              unsigned int & n_col)
{
  // get first and one-past-last iterator for row
  SparseMatrix<double>::const_iterator matrix_iterator = matrix.begin(i);
  SparseMatrix<double>::const_iterator matrix_iterator_end = matrix.end(i);

  // compute number of entries in row and then allocate memory
  n_col = matrix_iterator_end - matrix_iterator;
  row_values.resize(n_col);
  row_indices.resize(n_col);

  // loop over columns in row
  for (unsigned int k = 0; matrix_iterator != matrix_iterator_end;
       ++matrix_iterator, ++k)
  {
    row_values[k] = matrix_iterator->value();
    row_indices[k] = matrix_iterator->column();
  }
}

/**
 * \brief Checks that the imposed bounds are satisfied for a time step.
 */
/*
template <int dim>
bool FCT<dim>::check_max_principle(
  const Vector<double> & new_solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const double & dt)
{
  // machine precision for floating point comparisons
  const double machine_tolerance = 1.0e-12;

  // now set new precision
  std::cout.precision(15);

  // check that each dof value is bounded by its neighbors
  bool local_max_principle_satisfied = true;
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // check if dof is a Dirichlet node or not - if it is, don't check max
    // principle
    if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) ==
        dirichlet_nodes.end())
    {
      double value_i = new_solution(i);
      // check lower bound
      if (value_i < solution_min(i) - machine_tolerance)
      {
        local_max_principle_satisfied = false;
        // determine which condition was violated
        std::cout << "      Max principle lower bound violated with dof " << i
                  << " of low-order solution: " << std::scientific << value_i
                  << " < " << solution_min(i) << std::endl;
        debug_max_principle_low_order(i, low_order_ss_matrix, dt);
      }
      // check upper bound
      if (value_i > solution_max(i) + machine_tolerance)
      {
        local_max_principle_satisfied = false;
        // determine which condition was violated
        std::cout << "      Max principle upper bound violated with dof " << i
                  << " of low-order solution: " << std::scientific << value_i
                  << " > " << solution_max(i) << std::endl;
        debug_max_principle_low_order(i, low_order_ss_matrix, dt);
      }
    }
  }

  // restore default precision and format
  std::cout.unsetf(std::ios_base::floatfield);

  // exit if max principle was violated
  if (not local_max_principle_satisfied)
  {
    std::cout << "Program terminated due to max-principle violation."
              << std::endl;
    std::exit(0);
  }

  // return boolean for satisfaction of DMP
  return local_max_principle_satisfied;
}
*/

/** \brief Debugging function used for determining why the maximum
 *         principle is failing for the low-order method;
 *         examines the conditions that are required to be met for
 *         the principle to be satisfied.
 */
/*
template <int dim>
void FCT<dim>::debug_max_principle_low_order(
  const unsigned int & i,
  const SparseMatrix<double> & low_order_ss_matrix,
  const double & dt)
{
  // flag to determine if no conditions were violated
  bool condition_violated = false;

  // get nonzero entries of row i of A
  std::vector<double> row_values;
  std::vector<unsigned int> row_indices;
  unsigned int n_col;
  get_matrix_row(low_order_ss_matrix, i, row_values, row_indices, n_col);

  // diagonal entry
  double Aii = low_order_ss_matrix(i, i);

  // check that system matrix diagonal elements are non-negative
  if (Aii < 0.0)
  {
    std::cout << "         DEBUG: diagonal element is negative: " << Aii
              << std::endl;
    condition_violated = true;
  }

  // loop over nonzero entries in row i to compute off-diagonal sum
  // for diagonal dominance check and also check that each off-diagonal
  // element is non-positive
  double off_diagonal_sum = 0.0;
  for (unsigned int k = 0; k < n_col; ++k)
  {
    unsigned int j = row_indices[k];
    double Aij = row_values[k];

    if (j != i)
    {
      // add to off-diagonal sum
      off_diagonal_sum += std::abs(Aij);
      // check that system matrix off-diagonal elements are non-positive
      if (Aij > 0.0)
      {
        std::cout << "         DEBUG: off-diagonal element (" << i << "," << j
                  << ") is positive: " << Aij << std::endl;
        condition_violated = true;
      }
    }
  }

  // check that system matrix is diagonally dominant
  if (Aii < off_diagonal_sum)
  {
    std::cout << "         DEBUG: row is not diagonally dominant: Aii = " << Aii
              << ", off-diagonal sum = " << off_diagonal_sum << std::endl;
    condition_violated = true;
  }

  // check that CFL condition is satisfied if transient problem
  // if (!parameters.is_steady_state)
  //{
  double cfl = dt / (*lumped_mass_matrix)(i, i) * low_order_ss_matrix(i, i);
  if (cfl > 1.0)
  {
    std::cout << "         DEBUG: row does not satisfy CFL condition: CFL = "
              << cfl << std::endl;
    condition_violated = true;
  }
  //}

  // report if no conditions were violated
  if (not condition_violated)
    std::cout << "         DEBUG: No checks returned flags; deeper debugging is "
                 "necessary."
              << std::endl;
}
*/

/**
 * \brief Check to see if the DMP was satisfied at all time steps.
 *
 * \return flag telling whether imposed bounds were satisfied at all steps
 */
template <int dim>
bool FCT<dim>::check_DMP_satisfied()
{
  return DMP_satisfied_at_all_steps;
}

/**
 * \brief Outputs bounds to files.
 */
/*
template<int dim>
void FCT<dim>::output_bounds(
  const PostProcessor<dim> &postprocessor,
  const std::string &description_string) const
{
  // create the strings for the output files
  std::stringstream DMP_min_ss;
  std::stringstream DMP_max_ss;
  DMP_min_ss << "DMPmin_" << description_string;
  DMP_max_ss << "DMPmax_" << description_string;
  std::string DMP_min_string = DMP_min_ss.str();
  std::string DMP_max_string = DMP_max_ss.str();

  // use a method from the post-processor class to output the bounds
  postprocessor.output_solution(solution_min, *dof_handler, DMP_min_string);
  postprocessor.output_solution(solution_max, *dof_handler, DMP_max_string);
}
*/

/**
 * \brief Synchronizes limiting coefficients for a support point pair by taking
 *        the minimum limiting coefficient of all solution components.
 */
template <int dim>
void FCT<dim>::synchronize_min()
{
  // loop over support points
  for (unsigned int k_i = 0; k_i < n_dofs_scalar; ++k_i)
  {
    SparsityPattern::iterator it = scalar_sparsity_pattern.begin(k_i);
    SparsityPattern::iterator it_end = scalar_sparsity_pattern.end(k_i);
    for (; it != it_end; ++it)
    {
      // support point index
      const unsigned int k_j = it->column();

      // compute minimum limiting coefficient at this support point pair
      double min_ij = 1.0e15;
      for (unsigned int m = 0; m < n_components; ++m)
      {
        const unsigned int i = k_i * n_components + m;
        const unsigned int j = k_j * n_components + m;
        min_ij = std::min(min_ij, limiter_matrix(i, j));
      }

      // set limiting coefficients for this support point pair
      for (unsigned int m = 0; m < n_components; ++m)
      {
        const unsigned int i = k_i * n_components + m;
        const unsigned int j = k_j * n_components + m;
        limiter_matrix.set(i, j, min_ij);
      }
    }
  }
}

/**
 * \brief Computes the limited flux correction sum to ensure it is zero.
 *
 * \param[in] flux_vector limited correction flux sum vector
 */
template <int dim>
void FCT<dim>::check_limited_flux_correction_sum(
  const Vector<double> & flux_vector)
{
  // compute limited flux sum
  double sum = 0.0;
  for (unsigned int i = 0; i < n_dofs; ++i)
    sum += flux_vector[i];

  std::cout << std::scientific << std::setprecision(4)
            << "Limited flux sum = " << sum << std::endl;
}

/**
 * \brief Outputs the matrix of limiting coefficients.
 */
template <int dim>
void FCT<dim>::output_limiter_matrix() const
{
  // save limiting coefficients
  std::ofstream limiter_out("output/limiter.txt");
  limiter_matrix.print_formatted(
    limiter_out, 10, true, 0, "0", 1);
  limiter_out.close();
}
