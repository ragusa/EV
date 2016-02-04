/**
 * \brief Constructor for transient schemes.
 */
template <int dim>
FCT<dim>::FCT(const DoFHandler<dim> & dof_handler,
              Triangulation<dim> & triangulation,
              const SparseMatrix<double> & lumped_mass_matrix,
              const SparseMatrix<double> & consistent_mass_matrix,
              const LinearSolver<dim> & linear_solver,
              const SparsityPattern & sparsity_pattern,
              const std::vector<unsigned int> & dirichlet_nodes,
              const unsigned int & n_dofs,
              const unsigned int & dofs_per_cell,
              const bool & do_not_limit,
              const double & theta)
  : dof_handler(&dof_handler),
    triangulation(&triangulation),
    lumped_mass_matrix(&lumped_mass_matrix),
    consistent_mass_matrix(&consistent_mass_matrix),
    linear_solver(linear_solver),
    dirichlet_nodes(dirichlet_nodes),
    n_dofs(n_dofs),
    dofs_per_cell(dofs_per_cell),
    do_not_limit(do_not_limit),
    DMP_satisfied(true),
    theta(theta)
{
  // allocate memory for vectors
  tmp_vector.reinit(n_dofs);
  system_rhs.reinit(n_dofs);
  solution_min.reinit(n_dofs);
  solution_max.reinit(n_dofs);
  solution_min_new.reinit(n_dofs);
  solution_max_new.reinit(n_dofs);
  flux_correction_vector.reinit(n_dofs);
  Q_minus.reinit(n_dofs);
  Q_plus.reinit(n_dofs);
  R_minus.reinit(n_dofs);
  R_plus.reinit(n_dofs);

  // initialize sparse matrices
  system_matrix.reinit(sparsity_pattern);
  flux_correction_matrix.reinit(sparsity_pattern);
}

/**
 * \brief Constructor for steady-state schemes.
 */
template <int dim>
FCT<dim>::FCT(const DoFHandler<dim> & dof_handler,
              Triangulation<dim> & triangulation,
              const LinearSolver<dim> & linear_solver,
              const SparsityPattern & sparsity_pattern,
              const std::vector<unsigned int> & dirichlet_nodes,
              const unsigned int & n_dofs,
              const unsigned int & dofs_per_cell,
              const bool & do_not_limit)
  : dof_handler(&dof_handler),
    triangulation(&triangulation),
    linear_solver(linear_solver),
    dirichlet_nodes(dirichlet_nodes),
    n_dofs(n_dofs),
    dofs_per_cell(dofs_per_cell),
    do_not_limit(do_not_limit),
    DMP_satisfied(true),
    theta(0.0)
{
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
  system_matrix.reinit(sparsity_pattern);
  flux_correction_matrix.reinit(sparsity_pattern);
}

template <int dim>
void FCT<dim>::solve_FCT_system_fe(
  Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs,
  const double & dt,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const SparseMatrix<double> & high_order_diffusion_matrix)
{
  // compute flux corrections
  compute_flux_corrections_fe(new_solution,
                              old_solution,
                              dt,
                              low_order_diffusion_matrix,
                              high_order_diffusion_matrix);

  // compute max principle min and max values
  compute_bounds_fe(old_solution, low_order_ss_matrix, ss_rhs, dt);

  // compute limited flux bounds Q- and Q+
  compute_limited_flux_bounds_fe(old_solution, low_order_ss_matrix, ss_rhs, dt);

  // compute limited flux correction sum
  compute_limited_fluxes();

  // form transient rhs: system_rhs = M*u_old + dt*(ss_rhs - (A+D)*u_old + f)
  system_rhs = 0;
  system_rhs.add(dt, ss_rhs); //       now, system_rhs = dt*(ss_rhs)
  lumped_mass_matrix->vmult(tmp_vector, old_solution);
  system_rhs.add(1.0, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs)
  low_order_ss_matrix.vmult(tmp_vector, old_solution);
  system_rhs.add(
    -dt, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs - (A+D)*u_old)
  system_rhs.add(dt, flux_correction_vector); // now, system_rhs is complete

  // solve the linear system M*u_new = system_rhs
  system_matrix.copy_from(*lumped_mass_matrix);
  linear_solver.solve(system_matrix, new_solution, system_rhs, true);

  // check that local discrete maximum principle is satisfied at all time steps
  bool DMP_satisfied_this_step =
    check_max_principle(new_solution, low_order_ss_matrix, dt);
  DMP_satisfied = DMP_satisfied and DMP_satisfied_this_step;
}

/**
 * \brief Computes solution bounds for forward Euler time discretization.
 */
template <int dim>
void FCT<dim>::compute_bounds_fe(const Vector<double> & old_solution,
                                 const SparseMatrix<double> & low_order_ss_matrix,
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

    // find min and max values on cell
    double max_cell = old_solution(local_dof_indices[0]);
    double min_cell = old_solution(local_dof_indices[0]);
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      double value_j = old_solution(local_dof_indices[j]);
      max_cell = std::max(max_cell, value_j);
      min_cell = std::min(min_cell, value_j);
    }

    // update the max and min values of neighborhood of each dof
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      unsigned int i = local_dof_indices[j]; // global index
      solution_max(i) = std::max(solution_max(i), max_cell);
      solution_min(i) = std::min(solution_min(i), min_cell);
    }
  }

  // At this point, the min/max values of the old solution in the support
  // of test function i are stored in solution_min(i) and solution_max(i).
  // Now these values are multiplied by (1-dt/m(i))*sum_j(A(i,j)) and
  // added to dt/m(i)*b(i).

  // compute the upper and lower bounds for the maximum principle
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // compute sum of A_ij over row i
    // get nonzero entries of row i of A
    std::vector<double> row_values;
    std::vector<unsigned int> row_indices;
    unsigned int n_col;
    get_matrix_row(low_order_ss_matrix, i, row_values, row_indices, n_col);
    // add nonzero entries to get the row sum
    double row_sum = 0.0;
    for (unsigned int k = 0; k < n_col; ++k)
      row_sum += row_values[k];

    // compute the max and min values for the maximum principle
    solution_max(i) =
      solution_max(i) * (1.0 - dt / (*lumped_mass_matrix)(i, i) * row_sum) +
      dt / (*lumped_mass_matrix)(i, i) * ss_rhs(i);
    solution_min(i) =
      solution_min(i) * (1.0 - dt / (*lumped_mass_matrix)(i, i) * row_sum) +
      dt / (*lumped_mass_matrix)(i, i) * ss_rhs(i);
  }
}

/**
 * \brief Computes solution bounds for the transient case using a theta scheme.
 */
template <int dim>
void FCT<dim>::compute_bounds_theta(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs_new,
  const Vector<double> & ss_rhs_old,
  const double & dt)
{
  // references
  Vector<double> & solution_max_old = solution_max;
  Vector<double> & solution_min_old = solution_min;

  // initialize min and max of solution
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    solution_max_new(i) = new_solution(i);
    solution_min_new(i) = new_solution(i);
    solution_max_old(i) = old_solution(i);
    solution_min_old(i) = old_solution(i);
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

    // find min and max values on cell
    double max_cell_new = new_solution(local_dof_indices[0]);
    double min_cell_new = new_solution(local_dof_indices[0]);
    double max_cell_old = old_solution(local_dof_indices[0]);
    double min_cell_old = old_solution(local_dof_indices[0]);
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      double value_j_new = new_solution(local_dof_indices[j]);
      max_cell_new = std::max(max_cell_new, value_j_new);
      min_cell_new = std::min(min_cell_new, value_j_new);

      double value_j_old = old_solution(local_dof_indices[j]);
      max_cell_old = std::max(max_cell_old, value_j_old);
      min_cell_old = std::min(min_cell_old, value_j_old);
    }

    // update the max and min values of neighborhood of each dof
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      unsigned int i = local_dof_indices[j]; // global index
      solution_max_new(i) = std::max(solution_max_new(i), max_cell_new);
      solution_min_new(i) = std::min(solution_min_new(i), min_cell_new);
      solution_max_old(i) = std::max(solution_max_old(i), max_cell_old);
      solution_min_old(i) = std::min(solution_min_old(i), min_cell_old);
    }
  }

  // compute the upper and lower bounds for the maximum principle
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // compute sum of A_ij over row i
    // get nonzero entries of row i of A
    std::vector<double> row_values;
    std::vector<unsigned int> row_indices;
    unsigned int n_col;
    get_matrix_row(low_order_ss_matrix, i, row_values, row_indices, n_col);
    // add nonzero entries to get the row sum
    double off_diagonal_row_sum = 0.0;
    double diagonal_term = 0.0;
    for (unsigned int k = 0; k < n_col; ++k)
      if (row_indices[k] == i)
        diagonal_term = row_values[k];
      else
        off_diagonal_row_sum += row_values[k];
    // compute full row sum
    double row_sum = off_diagonal_row_sum + diagonal_term;

    // lumped mass matrix entry
    const double m_i = (*lumped_mass_matrix)(i, i);

    // compute the max and min values for the maximum principle
    solution_max(i) =
      ((1.0 - (1.0 - theta) * dt / m_i * row_sum) * solution_max_old(i) -
       theta * dt / m_i * off_diagonal_row_sum * solution_max_new(i) +
       dt / m_i * ((1 - theta) * ss_rhs_old(i) + theta * ss_rhs_new(i))) /
      (1.0 + theta * dt / m_i * diagonal_term);
    solution_min(i) =
      ((1.0 - (1.0 - theta) * dt / m_i * row_sum) * solution_min_old(i) -
       theta * dt / m_i * off_diagonal_row_sum * solution_min_new(i) +
       dt / m_i * ((1 - theta) * ss_rhs_old(i) + theta * ss_rhs_new(i))) /
      (1.0 + theta * dt / m_i * diagonal_term);
  }
}

/**
 * \brief Computes solution bounds for the steady-state case.
 */
template <int dim>
void FCT<dim>::compute_bounds_ss(const Vector<double> & solution,
                                 const SparseMatrix<double> & low_order_ss_matrix,
                                 const Vector<double> & ss_rhs)
{
  // initialize min and max values
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    solution_max(i) = solution(i);
    solution_min(i) = solution(i);
  }

  // loop over cells to get min and max solution values
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler->begin_active(),
                                                 endc = dof_handler->end();
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // find min and max values on cell
    double max_cell = solution(local_dof_indices[0]);
    double min_cell = solution(local_dof_indices[0]);
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      double value_j = solution(local_dof_indices[j]);
      max_cell = std::max(max_cell, value_j);
      min_cell = std::min(min_cell, value_j);
    }

    // update the max and min values of neighborhood of each dof
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      unsigned int i = local_dof_indices[j]; // global index
      solution_max(i) = std::max(solution_max(i), max_cell);
      solution_min(i) = std::min(solution_min(i), min_cell);
    }
  }

  // compute the upper and lower bounds for the maximum principle
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // get nonzero entries of row i of A
    std::vector<double> row_values;
    std::vector<unsigned int> row_indices;
    unsigned int n_col;
    get_matrix_row(low_order_ss_matrix, i, row_values, row_indices, n_col);

    // add nonzero entries to get the row sum
    double diagonal_term = 0.0;
    double off_diagonal_sum = 0.0;
    for (unsigned int k = 0; k < n_col; ++k)
      if (row_indices[k] == i)
        diagonal_term = row_values[k];
      else
        off_diagonal_sum += row_values[k];

    // compute the max and min values for the maximum principle
    solution_max(i) = -off_diagonal_sum / diagonal_term * solution_max(i) +
      ss_rhs(i) / diagonal_term;
    solution_min(i) = -off_diagonal_sum / diagonal_term * solution_min(i) +
      ss_rhs(i) / diagonal_term;
  }
}

/**
 * \brief Assembles the flux correction matrix for forward Euler
 *        time discretization.
 *
 * \param[in] dt current time step size
 */
template <int dim>
void FCT<dim>::compute_flux_corrections_fe(
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
 * \brief Assembles the flux correction matrix for transient schemes using a
 *        theta time discretization scheme.
 *
 * \param[in] dt current time step size
 */
template <int dim>
void FCT<dim>::compute_flux_corrections_theta(
  const Vector<double> & high_order_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const SparseMatrix<double> & high_order_diffusion_matrix_new,
  const SparseMatrix<double> & high_order_diffusion_matrix_old)
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
  SparseMatrix<double>::const_iterator it_high_new =
    high_order_diffusion_matrix_new.begin();
  SparseMatrix<double>::const_iterator it_high_old =
    high_order_diffusion_matrix_old.begin();
  SparseMatrix<double>::iterator it_flux = flux_correction_matrix.begin();
  SparseMatrix<double>::iterator it_end = flux_correction_matrix.end();

  for (; it_flux != it_end;
       ++it_flux, ++it_mass, ++it_low, ++it_high_new, ++it_high_old)
  {
    // get row and column indices
    const unsigned int i = it_flux->row();
    const unsigned int j = it_flux->column();

    // get values
    const double Mij = it_mass->value();
    const double DLij = it_low->value();
    const double DHij_new = it_high_new->value();
    const double DHij_old = it_high_old->value();

    // compute flux correction entry
    const double Fij = -Mij * (dUdt(j) - dUdt(i)) +
      (1.0 - theta) * (DLij - DHij_old) * (old_solution(j) - old_solution(i)) +
      theta * (DLij - DHij_new) *
        (high_order_solution(j) - high_order_solution(i));

    // store value
    it_flux->value() = Fij;
  }
}

/**
 * \brief Assembles the flux correction matrix for steady-state schemes.
 */
template <int dim>
void FCT<dim>::compute_flux_corrections_ss(
  const Vector<double> & high_order_solution,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const SparseMatrix<double> & high_order_diffusion_matrix)
{
  // reset flux correction matrix to zero
  flux_correction_matrix = 0;

  // iterate over sparse matrix entries
  SparseMatrix<double>::const_iterator it_low =
    low_order_diffusion_matrix.begin();
  SparseMatrix<double>::const_iterator it_high =
    high_order_diffusion_matrix.begin();
  SparseMatrix<double>::iterator it_flux = flux_correction_matrix.begin();
  SparseMatrix<double>::iterator it_end = flux_correction_matrix.end();

  for (; it_flux != it_end; ++it_flux, ++it_low, ++it_high)
  {
    // get row and column indices
    const unsigned int i = it_flux->row();
    const unsigned int j = it_flux->column();

    // get values
    const double DLij = it_low->value();
    const double DHij = it_high->value();

    // compute flux correction entry
    const double Fij =
      (DLij - DHij) * (high_order_solution(j) - high_order_solution(i));

    // store value
    it_flux->value() = Fij;
  }
}

/**
 * \brief Computes the limited flux bounds \f$\mathbf{Q}^+\f$ and
 *        \f$\mathbf{Q}^-\f$ for forward Euler time discretization.
 */
template <int dim>
void FCT<dim>::compute_limited_flux_bounds_fe(
  const Vector<double> & old_solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs,
  const double & dt)
{
  // compute Q+
  Q_plus = 0;
  lumped_mass_matrix->vmult(tmp_vector, solution_max);
  Q_plus.add(1.0 / dt, tmp_vector);
  lumped_mass_matrix->vmult(tmp_vector, old_solution);
  Q_plus.add(-1.0 / dt, tmp_vector);
  low_order_ss_matrix.vmult(tmp_vector, old_solution);
  Q_plus.add(1.0, tmp_vector);
  Q_plus.add(-1.0, ss_rhs);

  // compute Q-
  Q_minus = 0;
  lumped_mass_matrix->vmult(tmp_vector, solution_min);
  Q_minus.add(1.0 / dt, tmp_vector);
  lumped_mass_matrix->vmult(tmp_vector, old_solution);
  Q_minus.add(-1.0 / dt, tmp_vector);
  low_order_ss_matrix.vmult(tmp_vector, old_solution);
  Q_minus.add(1.0, tmp_vector);
  Q_minus.add(-1.0, ss_rhs);
}

/**
 * \brief Computes the limited flux bounds \f$\mathbf{Q}^+\f$ and
 *        \f$\mathbf{Q}^-\f$ for a theta time discretization.
 */
template <int dim>
void FCT<dim>::compute_limited_flux_bounds_theta(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs_new,
  const Vector<double> & ss_rhs_old,
  const double & dt)
{
  // compute Q+
  Q_plus = 0;
  lumped_mass_matrix->vmult(tmp_vector, solution_max);
  Q_plus.add(1.0 / dt, tmp_vector);
  lumped_mass_matrix->vmult(tmp_vector, old_solution);
  Q_plus.add(-1.0 / dt, tmp_vector);
  low_order_ss_matrix.vmult(tmp_vector, old_solution);
  Q_plus.add(1.0 - theta, tmp_vector);
  low_order_ss_matrix.vmult(tmp_vector, new_solution);
  Q_plus.add(theta, tmp_vector);
  Q_plus.add(-(1.0 - theta), ss_rhs_old);
  Q_plus.add(-theta, ss_rhs_new);

  // compute Q-
  Q_minus = 0;
  lumped_mass_matrix->vmult(tmp_vector, solution_min);
  Q_minus.add(1.0 / dt, tmp_vector);
  lumped_mass_matrix->vmult(tmp_vector, old_solution);
  Q_minus.add(-1.0 / dt, tmp_vector);
  low_order_ss_matrix.vmult(tmp_vector, old_solution);
  Q_minus.add(1.0 - theta, tmp_vector);
  low_order_ss_matrix.vmult(tmp_vector, new_solution);
  Q_minus.add(theta, tmp_vector);
  Q_minus.add(-(1.0 - theta), ss_rhs_old);
  Q_minus.add(-theta, ss_rhs_new);
}

/**
 * \brief Computes the limited flux bounds \f$\mathbf{Q}^+\f$ and
 *        \f$\mathbf{Q}^-\f$ for the steady-state case.
 *
 * \f[
 *   Q_i^\pm \equiv A_{i,i}^L W_i^\pm + \sum\limits_{j\ne i}A_{i,j}^L U_j - b_i
 * \f]
 *
 * \param[in] solution solution iterate \f$\mathbf{U}^{(\ell)}\f$
 * \param[in] low_order_ss_matrix low-order steady-state matrix \f$\mathbf{A}^L\f$
 * \param[in] ss_rhs steady-state right-hand-side vector \f$\mathbf{b}\f$
 */
template <int dim>
void FCT<dim>::compute_limited_flux_bounds_ss(
  const Vector<double> & solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs)
{
  // compute Q- and Q+
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // get nonzero entries in row i of low-order steady-state matrix
    std::vector<double> row_values;
    std::vector<unsigned int> row_indices;
    unsigned int n_col;
    get_matrix_row(low_order_ss_matrix, i, row_values, row_indices, n_col);

    // compute Q-(i) and Q+(i)
    Q_minus[i] = -ss_rhs[i];
    Q_plus[i] = -ss_rhs[i];
    for (unsigned int k = 0; k < n_col; ++k)
    {
      if (row_indices[k] == i) // diagonal element of A^L
      {
        Q_minus[i] += row_values[k] * solution_min[i];
        Q_plus[i] += row_values[k] * solution_max[i];
      }
      else // off-diagonal element of A^L
      {
        const unsigned int j = row_indices[k];
        Q_minus[i] += row_values[k] * solution[j];
        Q_plus[i] += row_values[k] * solution[j];
      }
    }
  }
}

/**
 * \brief Computes the limited flux vectors for any FCT scheme.
 */
template <int dim>
void FCT<dim>::compute_limited_fluxes()
{
  // reset flux correction vector
  flux_correction_vector = 0;

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

  // if user chose not to limit, set R+ and R- = 1 so that limiting
  // coefficients will equal 1 and thus no limiting will occur
  if (do_not_limit)
  {
    for (unsigned int i = 0; i < n_dofs; ++i)
    {
      R_minus(i) = 1.0;
      R_plus(i) = 1.0;
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

      // add Lij*Fij to flux correction sum
      flux_correction_vector(i) += Lij * Fij;
    }
  }
}

/** \brief Gets the values and indices of nonzero elements in a row of a sparse
 * matrix.
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
 * \brief Check that the DMP is satisfied for a time step.
 */
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

/** \brief Debugging function used for determining why the maximum
 *         principle is failing for the low-order method;
 *         examines the conditions that are required to be met for
 *         the principle to be satisfied.
 */
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

/**
 * \brief Check to see if the DMP was satisfied at all time steps.
 */
template <int dim>
bool FCT<dim>::check_DMP_satisfied()
{
  return DMP_satisfied;
}

/** \brief Outputs bounds to files.
 */
template <int dim>
void FCT<dim>::output_bounds(const PostProcessor<dim> & postprocessor,
                             const std::string & description_string) const
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

/**
 * \brief Returns the limited flux correction vector
 *
 * \return the limited flux correction vector
 */
template <int dim>
Vector<double> FCT<dim>::get_flux_correction_vector() const
{
  return flux_correction_vector;
}
