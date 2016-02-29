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
              const bool & include_analytic_bounds_,
              const FESystem<dim> & fe_,
              const QGauss<dim> & cell_quadrature_,
              const FunctionParser<dim> & cross_section_function_,
              FunctionParser<dim> & source_function_,
              const bool & source_is_time_dependent_,
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
    include_analytic_bounds(include_analytic_bounds_),
    fe(&fe_),
    cell_quadrature(&cell_quadrature_),
    n_q_points_cell(cell_quadrature_.size()),
    cross_section_function(&cross_section_function_),
    source_function(&source_function_),
    source_is_time_dependent(source_is_time_dependent_),
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
  reaction_min.reinit(n_dofs);
  reaction_max.reinit(n_dofs);
  source_min.reinit(n_dofs);
  source_max.reinit(n_dofs);
  lower_bound_analytic.reinit(n_dofs);
  upper_bound_analytic.reinit(n_dofs);
  flux_correction_vector.reinit(n_dofs);
  Q_minus.reinit(n_dofs);
  Q_plus.reinit(n_dofs);
  R_minus.reinit(n_dofs);
  R_plus.reinit(n_dofs);

  // initialize sparse matrices
  system_matrix.reinit(sparsity_pattern);
  flux_correction_matrix.reinit(sparsity_pattern);
  limited_flux_matrix.reinit(sparsity_pattern);

  // compute min and max of reaction coefficients in neighborhood of each DoF
  compute_function_bounds(*cross_section_function, reaction_min, reaction_max);

  // compute min and max of source in neighborhood of each DoF
  if (!source_is_time_dependent)
    compute_function_bounds(*source_function, source_min, source_max);
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
              const bool & do_not_limit,
              const bool & include_analytic_bounds_,
              const FESystem<dim> & fe_,
              const QGauss<dim> & cell_quadrature_,
              const FunctionParser<dim> & cross_section_function_,
              FunctionParser<dim> & source_function_)
  : dof_handler(&dof_handler),
    triangulation(&triangulation),
    linear_solver(linear_solver),
    dirichlet_nodes(dirichlet_nodes),
    n_dofs(n_dofs),
    dofs_per_cell(dofs_per_cell),
    do_not_limit(do_not_limit),
    include_analytic_bounds(include_analytic_bounds_),
    fe(&fe_),
    cell_quadrature(&cell_quadrature_),
    n_q_points_cell(cell_quadrature_.size()),
    cross_section_function(&cross_section_function_),
    source_function(&source_function_),
    source_is_time_dependent(false),
    DMP_satisfied(true),
    theta(0.0)
{
  // allocate memory for vectors
  tmp_vector.reinit(n_dofs);
  system_rhs.reinit(n_dofs);
  solution_min.reinit(n_dofs);
  solution_max.reinit(n_dofs);
  reaction_min.reinit(n_dofs);
  reaction_max.reinit(n_dofs);
  source_min.reinit(n_dofs);
  source_max.reinit(n_dofs);
  lower_bound_analytic.reinit(n_dofs);
  upper_bound_analytic.reinit(n_dofs);
  flux_correction_vector.reinit(n_dofs);
  Q_minus.reinit(n_dofs);
  Q_plus.reinit(n_dofs);
  R_minus.reinit(n_dofs);
  R_plus.reinit(n_dofs);

  // initialize sparse matrices
  system_matrix.reinit(sparsity_pattern);
  flux_correction_matrix.reinit(sparsity_pattern);
  limited_flux_matrix.reinit(sparsity_pattern);

  // compute min and max of reaction coefficients in neighborhood of each DoF
  compute_function_bounds(*cross_section_function, reaction_min, reaction_max);

  // compute min and max of source in neighborhood of each DoF
  compute_function_bounds(*source_function, source_min, source_max);
}

template <int dim>
void FCT<dim>::solve_FCT_system_fe(
  Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const SparseMatrix<double> & low_order_ss_matrix,
  const Vector<double> & ss_rhs,
  const double & dt,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const SparseMatrix<double> & high_order_diffusion_matrix,
  const double & t_old)
{
  // compute flux corrections
  compute_flux_corrections_fe(new_solution,
                              old_solution,
                              dt,
                              low_order_diffusion_matrix,
                              high_order_diffusion_matrix);

  // compute max principle min and max values
  compute_bounds_fe(old_solution, low_order_ss_matrix, ss_rhs, dt, t_old);

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
  bool DMP_satisfied_this_step = check_fct_bounds(new_solution);
  DMP_satisfied = DMP_satisfied and DMP_satisfied_this_step;
}

/**
 * \brief Computes solution bounds for forward Euler time discretization.
 */
template <int dim>
void FCT<dim>::compute_bounds_fe(const Vector<double> & old_solution,
                                 const SparseMatrix<double> & low_order_ss_matrix,
                                 const Vector<double> & ss_rhs,
                                 const double & dt,
                                 const double & time_old)
{
  // compute min and max of solution in neighborhood of each DoF
  compute_min_max_solution(old_solution, solution_min, solution_max);

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

  // apply analytic bounds if specified
  if (include_analytic_bounds)
    compute_and_apply_analytic_bounds(old_solution, dt, time_old);
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
  const double & dt,
  const double & time_old)
{
  // references
  Vector<double> & solution_max_old = solution_max;
  Vector<double> & solution_min_old = solution_min;

  // compute min and max of solution in neighborhood of each DoF
  compute_min_max_solution(old_solution, solution_min_old, solution_max_old);
  compute_min_max_solution(new_solution, solution_min_new, solution_max_new);

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

  // apply analytic bounds if specified
  if (include_analytic_bounds)
    compute_and_apply_analytic_bounds(old_solution, dt, time_old);
}

/**
 * \brief Computes solution bounds for the steady-state case.
 */
template <int dim>
void FCT<dim>::compute_bounds_ss(const Vector<double> & solution,
                                 const SparseMatrix<double> & low_order_ss_matrix,
                                 const Vector<double> & ss_rhs)
{
  // compute min and max of solution in neighborhood of each DoF
  compute_min_max_solution(solution, solution_min, solution_max);

  // compute the upper and lower bounds for the maximum principle
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) ==
        dirichlet_nodes.end())
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
    else
    {
      solution_max(i) = solution(i);
      solution_min(i) = solution(i);
    }
  }

  // apply analytic bounds if specified
  if (include_analytic_bounds)
    compute_and_apply_analytic_bounds(solution, 0.0, 0.0);
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
  const Vector<double> & cumulative_antidiffusion,
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
  Q_plus.add(-1.0, cumulative_antidiffusion);

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
  Q_minus.add(-1.0, cumulative_antidiffusion);
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
  // compute the limited flux correction matrix: L * P
  compute_limited_flux_matrix();

  // compute the row-sum of the limited flux correction matrix: row_sum(L * P)
  compute_row_sums(limited_flux_matrix, flux_correction_vector);
}

/**
 * \brief Computes the limited flux matrix for any FCT scheme.
 */
template <int dim>
void FCT<dim>::compute_limited_flux_matrix()
{
  // reset limited flux matrix
  limited_flux_matrix = 0;

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

  // compute limited flux matrix
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // get first and one-past-last iterator for row
    SparseMatrix<double>::const_iterator it_flux =
      flux_correction_matrix.begin(i);
    SparseMatrix<double>::iterator it_limflux = limited_flux_matrix.begin(i);
    SparseMatrix<double>::const_iterator it_flux_end =
      flux_correction_matrix.end(i);

    // loop over columns in row
    for (; it_flux != it_flux_end; ++it_flux, ++it_limflux)
    {
      // get column index
      const unsigned int j = it_flux->column();

      // get flux value
      const double Pij = it_flux->value();

      // compute limiting coefficient
      double Lij;
      if (Pij >= 0.0)
        Lij = std::min(R_plus(i), R_minus(j));
      else
        Lij = std::min(R_minus(i), R_plus(j));

      // store limited flux
      it_limflux->value() = Lij * Pij;
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
bool FCT<dim>::check_fct_bounds(const Vector<double> & solution) const
{
  // machine precision for floating point comparisons
  const double machine_tolerance = 1.0e-15;

  // now set new precision
  std::cout.precision(15);

  // check that each dof value is bounded by its neighbors
  bool fct_bounds_satisfied = true;
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // check bounds if dof does not correspond to a Dirichlet node
    if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) ==
        dirichlet_nodes.end())
    {
      double value_i = solution(i);

      // check lower bound
      if (value_i < solution_min(i) - machine_tolerance)
      {
        fct_bounds_satisfied = false;

        std::cout << "FCT bounds violated by dof " << i << ": " << value_i
                  << " < " << solution_min(i) << std::endl;
      }
      // check upper bound
      if (value_i > solution_max(i) + machine_tolerance)
      {
        fct_bounds_satisfied = false;

        std::cout << "FCT bounds violated by dof " << i << ": " << value_i
                  << " > " << solution_max(i) << std::endl;
      }
    }
  }

  // restore default precision and format
  std::cout.unsetf(std::ios_base::floatfield);

  // return boolean for satisfaction of DMP
  return fct_bounds_satisfied;
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
void FCT<dim>::output_bounds(const PostProcessor<dim> & postprocessor) const
{
  // use a method from the post-processor class to output the bounds
  postprocessor.output_solution(solution_min, *dof_handler, "DMPmin");
  postprocessor.output_solution(solution_max, *dof_handler, "DMPmax");
}

/**
 * \brief Returns the limited flux correction vector
 *
 * \return the limited flux correction vector
 */
template <int dim>
Vector<double> FCT<dim>::get_limited_flux_vector() const
{
  return flux_correction_vector;
}

/**
 * \brief Computes the minimum and maximum of a time-dependent function
 *        in the neighborhood of each degree of freedom.
 *
 * \param[in] function     function
 * \param[in] time         time at which to evaluate function if time-dependent
 * \param[out] min_values  vector of minimum values in neighborhood of each DoF
 * \param[out] max_values  vector of maximum values in neighborhood of each DoF
 */
template <int dim>
void FCT<dim>::compute_function_bounds(FunctionParser<dim> & function,
                                       const double & time,
                                       Vector<double> & min_values,
                                       Vector<double> & max_values) const
{
  // set time for function
  function.set_time(time);

  // call time-independent version of function
  compute_function_bounds(function, min_values, max_values);
}

/**
 * \brief Computes the minimum and maximum of a function
 *        in the neighborhood of each degree of freedom.
 *
 * \param[in] function     function
 * \param[out] min_values  vector of minimum values in neighborhood of each DoF
 * \param[out] max_values  vector of maximum values in neighborhood of each DoF
 */
template <int dim>
void FCT<dim>::compute_function_bounds(const FunctionParser<dim> & function,
                                       Vector<double> & min_values,
                                       Vector<double> & max_values) const
{
  // initialize min and max values
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    min_values(i) = 1.0e15;
    max_values(i) = -1.0e15;
  }

  // FE values
  FEValues<dim> fe_values(*fe, *cell_quadrature, update_quadrature_points);

  // loop over cells
  Cell cell = dof_handler->begin_active(), endc = dof_handler->end();
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  for (; cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // get quadrature points on cell
    std::vector<Point<dim>> points(n_q_points_cell);
    points = fe_values.get_quadrature_points();

    // find min and max values on cell
    double min_cell = 1.0e15;
    double max_cell = -1.0e15;
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      // update min and max values on cell
      const double value_q = function.value(points[q]);
      min_cell = std::min(min_cell, value_q);
      max_cell = std::max(max_cell, value_q);
    }

    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // update min and max values in neighborhood of each dof
    for (unsigned int j = 0; j < dofs_per_cell; ++j)
    {
      unsigned int i = local_dof_indices[j]; // global index
      min_values(i) = std::min(min_values(i), min_cell);
      max_values(i) = std::max(max_values(i), max_cell);
    }
  }
}

/**
 * \brief Computes the minimum and maximum of the solution
 *        in the neighborhood of each degree of freedom.
 *
 * \param[in] solution       solution at which to evaluate min/max
 * \param[out] min_solution  min solution in neighborhood of each DoF
 * \param[out] max_solution  max solution in neighborhood of each DoF
 */
template <int dim>
void FCT<dim>::compute_min_max_solution(const Vector<double> & solution,
                                        Vector<double> & min_solution,
                                        Vector<double> & max_solution) const
{
  // initialize min and max
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    min_solution(i) = solution(i);
    max_solution(i) = solution(i);
  }

  // loop over cells
  Cell cell = dof_handler->begin_active(), endc = dof_handler->end();
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
      min_solution(i) = std::min(min_solution(i), min_cell);
      max_solution(i) = std::max(max_solution(i), max_cell);
    }
  }
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
void FCT<dim>::compute_and_apply_analytic_bounds(const Vector<double> & solution,
                                                 const double & dt,
                                                 const double & time_old)
{
  // compute min and max of source in neighborhood of each DoF if source
  // is time-dependent
  if (source_is_time_dependent)
    compute_function_bounds(*source_function, time_old, source_min, source_max);
  else
    compute_function_bounds(*source_function, source_min, source_max);

  // compute min and max of solution in neighborhood of each DoF
  compute_min_max_solution(solution, lower_bound_analytic, upper_bound_analytic);

  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // branch on minimum reaction coefficient to avoid division by zero
    double min_source_term;
    if (std::abs(reaction_max[i]) < 1.0e-15) // equal to zero within precision
      min_source_term = dt * source_min[i];
    else
      min_source_term =
        source_min[i] / reaction_max[i] * (1.0 - std::exp(-dt * reaction_max[i]));

    // branch on maximum reaction coefficient to avoid division by zero
    double max_source_term;
    if (std::abs(reaction_min[i]) < 1.0e-15) // equal to zero within precision
      max_source_term = dt * source_max[i];
    else
      max_source_term =
        source_max[i] / reaction_min[i] * (1.0 - std::exp(-dt * reaction_min[i]));

    // compute bounds
    lower_bound_analytic[i] =
      lower_bound_analytic[i] * std::exp(-dt * reaction_max[i]) + min_source_term;
    upper_bound_analytic[i] =
      upper_bound_analytic[i] * std::exp(-dt * reaction_min[i]) + max_source_term;

    // apply bounds
    solution_min[i] = std::min(solution_min[i], lower_bound_analytic[i]);
    solution_max[i] = std::max(solution_max[i], upper_bound_analytic[i]);
  }
}

/**
 * \brief Computes the row-sum vector of a matrix.
 *
 * \param[in]  matrix   matrix
 * \param[out] row_sum  vector of row-sums
 */
template <int dim>
void FCT<dim>::compute_row_sums(const SparseMatrix<double> & matrix,
                                Vector<double> & row_sums) const
{
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    // get first and one-past-last iterator for row
    SparseMatrix<double>::const_iterator matrix_iterator = matrix.begin(i);
    SparseMatrix<double>::const_iterator matrix_iterator_end = matrix.end(i);

    // compute row sum
    row_sums(i) = 0.0;
    for (; matrix_iterator != matrix_iterator_end; ++matrix_iterator)
      row_sums(i) += matrix_iterator->value();
  }
}

/**
 * \brief Subtracts the limited flux correction matrix from the remaining flux
 *        correction matrix.
 *
 * This function is to be used with the cumulative antidiffusive flux algorithm
 * for implicit schemes.
 */
template <int dim>
void FCT<dim>::subtract_limited_flux_correction_matrix()
{
  flux_correction_matrix.add(-1.0, limited_flux_matrix);
}
