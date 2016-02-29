/**
 * \file FCT.cc
 * \brief Provides the function definitions for the FCT class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] parameters_  parameters
 * \param[in] dof_handler_  DoF handler
 * \param[in] triangulation_  triangulation
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 * \param[in] consistent_mass_matrix_  consistent mass matrix \f$\mathbf{M}^C\f$
 * \param[in] star_state_  star state object for computing
 *                         \f$\mathbf{U}^*_{i,j}\f$
 * \param[in] linear_solver_  linear solver
 * \param[in] sparsity_pattern_  sparsity pattern
 * \param[in] dirichlet_nodes_  vector of Dirichlet nodes
 * \param[in] n_components_  number of components
 * \param[in] dofs_per_cell_  number of degrees of freedom per cell
 * \param[in] component_names_  list of component names
 * \param[in] use_star_states_in_fct_bounds_  flag to signal that star states are
 *            to be used in FCT bounds
 */
template <int dim>
FCT<dim>::FCT(const RunParameters<dim> & parameters_,
              const DoFHandler<dim> & dof_handler_,
              const Triangulation<dim> & triangulation_,
              const SparseMatrix<double> & lumped_mass_matrix_,
              const SparseMatrix<double> & consistent_mass_matrix_,
              const std::shared_ptr<StarState<dim>> & star_state_,
              const LinearSolver<dim> & linear_solver_,
              const SparsityPattern & sparsity_pattern_,
              const std::vector<unsigned int> & dirichlet_nodes_,
              const unsigned int & n_components_,
              const unsigned int & dofs_per_cell_,
              const std::vector<std::string> & component_names_,
              const bool & use_star_states_in_fct_bounds_)
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
    fct_bounds_type(parameters_.fct_bounds_type),
    antidiffusion_type(parameters_.antidiffusion_type),
    synchronization_type(parameters_.fct_synchronization_type),
    fct_limitation_type(parameters_.fct_limitation_type),
    use_star_states_in_fct_bounds(use_star_states_in_fct_bounds_),
    fct_bounds_satisfied_at_all_steps(true),
    bounds_transient_file_index(1)
{
  // assert that DMP bounds are not used for the non-scalar case
  // if (fct_bounds_type == FCTBoundsType::dmp)
  //  Assert(n_components == 1, ExcNotImplemented());

  // reinitialize
  reinitialize(sparsity_pattern_);

  // create bounds component names
  lower_bound_component_names = component_names_;
  upper_bound_component_names = component_names_;
  for (unsigned int m = 0; m < n_components; ++m)
  {
    lower_bound_component_names[m] += "_FCTmin";
    upper_bound_component_names[m] += "_FCTmax";
  }
}

/**
 * \brief Reinitializes.
 *
 * \param[in] sparsity_pattern  sparsity pattern
 */
template <int dim>
void FCT<dim>::reinitialize(const SparsityPattern & sparsity_pattern)
{
  // distribute degrees of freedom in scalar degree of freedom handler
  dof_handler_scalar.clear();
  dof_handler_scalar.distribute_dofs(fe_scalar);

  // allocate memory for vectors
  old_solution_characteristic.reinit(n_dofs);
  tmp_vector.reinit(n_dofs);
  system_rhs.reinit(n_dofs);
  solution_min.reinit(n_dofs);
  solution_max.reinit(n_dofs);
  flux_correction_vector.reinit(n_dofs);
  Q_minus.reinit(n_dofs);
  Q_plus.reinit(n_dofs);
  R_minus.reinit(n_dofs);
  R_plus.reinit(n_dofs);

  // create scalar sparsity pattern
  DynamicSparsityPattern dsp(n_dofs_scalar);
  DoFTools::make_sparsity_pattern(dof_handler_scalar, dsp);
  scalar_sparsity_pattern.copy_from(dsp);

  // initialize sparse matrices
  system_matrix.reinit(sparsity_pattern);
  limiter_matrix.reinit(sparsity_pattern);
  flux_correction_matrix.reinit(sparsity_pattern);
  flux_correction_matrix_transformed.reinit(sparsity_pattern);

  // if needed, create lists of DoF indices; this is needed if limiting is
  // performed on a non-conservative set of variables such as characteristic
  create_dof_indices_lists();
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

  // branch on type of limitation: conservative variables or characteristic
  switch (fct_limitation_type)
  {
    case FCTLimitationType::conservative:
      // compute FCT solution bounds W- and W+
      compute_solution_bounds(old_solution, ss_reaction, ss_rhs, dt);

      // compute antidiffusion bounds Q- and Q+
      compute_antidiffusion_bounds(
        old_solution, ss_flux, ss_rhs, low_order_diffusion_matrix, dt);

      break;
    case FCTLimitationType::characteristic:
      // compute characteristic FCT solution bounds \hat{W}- and \hat{W}+
      compute_solution_bounds_characteristic(
        old_solution, ss_reaction, ss_rhs, dt);

      // transform antidiffusive fluxes: P_{i,j} -> \hat{P}_{i,j}
      transform_matrix(
        old_solution, flux_correction_matrix, flux_correction_matrix_transformed);

      // compute characteristic antidiffusion bounds \hat{Q}- and \hat{Q}+
      compute_antidiffusion_bounds_characteristic(old_solution_characteristic,
                                                  dt);

      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // compute limited flux correction sum and add it to rhs
  switch (antidiffusion_type)
  {
    case AntidiffusionType::limited: // normal, limited antidiffusion flux
      // compute limiter matrix
      compute_limiting_coefficients_zalesak();

      // synchronize limiting coefficients if requested
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

      // compute limited flux sums
      compute_limited_flux_correction_vector();

      break;
    case AntidiffusionType::full: // full antidiffusion flux
      compute_full_flux_correction_vector();
      break;
    case AntidiffusionType::none: // no antidiffusion flux
      flux_correction_vector = 0;
      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

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

  // check that local discrete maximum principle is satisfied at all time steps
  if (fct_limitation_type == FCTLimitationType::conservative &&
      antidiffusion_type != AntidiffusionType::full)
  {
    bool fct_bounds_satisfied_this_step =
      check_fct_bounds_satisfied(new_solution);

    fct_bounds_satisfied_at_all_steps =
      fct_bounds_satisfied_at_all_steps && fct_bounds_satisfied_this_step;
  }
}

/**
 * \brief Computes the minimum and maximum degree of freedom values in the
 *        support of each degree of freedom.
 *
 * \param[in] solution  solution vector \f$\mathbf{U}\f$
 * \param[out] min_values  minimum values of solution:
 *             \f$U_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} U_j\f$
 * \param[out] max_values  maximum values of solution:
 *             \f$U_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} U_j\f$
 */
template <int dim>
void FCT<dim>::compute_min_and_max_of_solution(const Vector<double> & solution,
                                               Vector<double> & min_values,
                                               Vector<double> & max_values) const
{
  // initialize solution min and max
  for (unsigned int i = 0; i < n_dofs; ++i)
  {
    min_values(i) = solution(i);
    max_values(i) = solution(i);
  }

  // loop over cells
  Cell cell = dof_handler->begin_active(), endc = dof_handler->end();
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // loop over solution components
    for (unsigned int m = 0; m < n_components; ++m)
    {
      // find min and max values on cell for this component
      double max_cell = solution(local_dof_indices[m]);
      double min_cell = solution(local_dof_indices[m]);
      for (unsigned int j = 0; j < dofs_per_cell_per_component; ++j)
      {
        // consider old states
        const double value_j = solution(local_dof_indices[j * n_components + m]);
        max_cell = std::max(max_cell, value_j);
        min_cell = std::min(min_cell, value_j);
      }

      // update the max and min values of neighborhood of each dof
      for (unsigned int j = 0; j < dofs_per_cell_per_component; ++j)
      {
        unsigned int i = local_dof_indices[j * n_components + m]; // global index
        max_values(i) = std::max(max_values(i), max_cell);
        min_values(i) = std::min(min_values(i), min_cell);
      }
    }
  }
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
 * where \f$U_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} U_j\f$ and
 * \f$U_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} U_j\f$.
 *
 * \param[in] old_solution  old solution vector \f$\mathbf{U}^n\f$
 * \param[in] ss_reaction  steady-state reaction vector \f$\mathbf{\sigma}\f$,
 *            where \f$
 *   \sigma_i = \int\limits_{S_i}\varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *     = \sum\limits_j\int\limits_{S_{i,j}}
 *     \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 */
template <int dim>
void FCT<dim>::compute_solution_bounds(const Vector<double> & old_solution,
                                       const Vector<double> & ss_reaction,
                                       const Vector<double> & ss_rhs,
                                       const double & dt)
{
  // compute minimum and maximum values of solution
  compute_min_and_max_of_solution(old_solution, solution_min, solution_max);

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

  // compute the upper and lower bounds for the FCT solution
  switch (fct_bounds_type)
  {
    case FCTBoundsType::dmp:
    {
      for (unsigned int i = 0; i < n_dofs; ++i)
      {
        solution_max(i) = solution_max(i) *
            (1.0 - dt / (*lumped_mass_matrix)(i, i) * ss_reaction(i)) +
          dt / (*lumped_mass_matrix)(i, i) * ss_rhs(i);

        solution_min(i) = solution_min(i) *
            (1.0 - dt / (*lumped_mass_matrix)(i, i) * ss_reaction(i)) +
          dt / (*lumped_mass_matrix)(i, i) * ss_rhs(i);
      }
      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }
}

/**
 * \brief Computes FCT solution bounds for use with the characteristic limiter,
 *        \f$\hat{W}_i^-\f$ and \f$\hat{W}_i^+\f$.
 *
 * The bounds are computed as
 * \f[
 *   \hat{W}_i^-(\mathrm{U}^n) = \hat{U}_{min,i}^n \,,
 * \f]
 * \f[
 *   \hat{W}_i^+(\mathrm{U}^n) = \hat{U}_{max,i}^n \,,
 * \f]
 * where \f$\hat{U}_{min,i}\equiv\min\limits_{j\in \mathcal{I}(S_i)} \hat{U}_j\f$
 * and \f$\hat{U}_{max,i}\equiv\max\limits_{j\in \mathcal{I}(S_i)} \hat{U}_j\f$.
 *
 * \warning These bounds currently assume
 * - no reaction term
 * - no source term
 *
 * \param[in] old_solution  old solution vector \f$\mathbf{U}^n\f$
 * \param[in] ss_reaction  steady-state reaction vector \f$\mathbf{\sigma}\f$,
 *            where \f$
 *   \sigma_i = \int\limits_{S_i}\varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *     = \sum\limits_j\int\limits_{S_{i,j}}
 *     \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 */
template <int dim>
void FCT<dim>::compute_solution_bounds_characteristic(
  const Vector<double> & old_solution,
  const Vector<double> &,
  const Vector<double> &,
  const double &)
{
  // transform old solution vector to characteristic variables
  transform_vector(old_solution, old_solution, old_solution_characteristic);

  // compute minimum and maximum values of solution
  compute_min_and_max_of_solution(
    old_solution_characteristic, solution_min, solution_max);
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
 * \brief Computes the antidiffusion bounds \f$\mathbf{Q}^\pm\f$.
 *
 * \param[in] old_solution  old solution \f$\mathbf{U}^n\f$
 * \param[in] ss_flux  steady-state inviscid flux vector
 *   (entries are \f$(\mathbf{A}\mathbf{U}^n)_i\f$ for scalar case,
 *   \f$\sum\limits_j\mathbf{c}_{i,j}\cdot\mathrm{F}^n_j\f$ for systems case)
 * \param[in] ss_rhs  steady-state right hand side vector \f$\mathbf{b}^n\f$
 * \param[in] low_order_diffusion_matrix  low-order diffusion matrix
 *   \f$\mathbf{D}^L\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 */
template <int dim>
void FCT<dim>::compute_antidiffusion_bounds(
  const Vector<double> & old_solution,
  const Vector<double> & ss_flux,
  const Vector<double> & ss_rhs,
  const SparseMatrix<double> & low_order_diffusion_matrix,
  const double & dt)
{
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
}

/**
 * \brief Computes the characteristic antidiffusion bounds
 *        \f$\hat{\mathbf{Q}}^\pm\f$.
 *
 * \param[in] old_solution  characteristic old solution \f$\hat{\mathbf{U}}^n\f$
 * \param[in] dt  time step size \f$\Delta t\f$
 */
template <int dim>
void FCT<dim>::compute_antidiffusion_bounds_characteristic(
  const Vector<double> & old_solution, const double & dt)
{
  // compute Q-
  Q_minus = 0;
  lumped_mass_matrix->vmult(tmp_vector, solution_min);
  Q_minus.add(1.0 / dt, tmp_vector);
  lumped_mass_matrix->vmult(tmp_vector, old_solution);
  Q_minus.add(-1.0 / dt, tmp_vector);

  // compute Q+
  Q_plus = 0;
  lumped_mass_matrix->vmult(tmp_vector, solution_max);
  Q_plus.add(1.0 / dt, tmp_vector);
  lumped_mass_matrix->vmult(tmp_vector, old_solution);
  Q_plus.add(-1.0 / dt, tmp_vector);
}

/**
 * \brief Computes the limiting coefficient matrix.
 *
 * This function computes the limiting coefficient matrix \f$\mathbf{L}\f$
 * and stores it in \link limiter_matrix \endlink.
 *
 * \pre This function assumes
 * - the FCT solution bounds \f$\mathbf{W}^\pm\f$
 *   have been computed and stored in \link solution_max \endlink and \link
 *   solution_min \endlink, respectively.
 * - the antidiffusive correction flux matrix \f$\mathbf{P}\f$ has been computed
 *   and stored in \link flux_correction_matrix \endlink.
 */
template <int dim>
void FCT<dim>::compute_limiting_coefficients_zalesak()
{
  // reset limiter matrix
  limiter_matrix = 0;

  // compute R- and R+
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

      // std::cout<<"L("<<i<<","<<j<<") = "<<Lij<<std::endl;

      // store entry in limiter matrix
      limiter_matrix.set(i, j, Lij);
    }
  }
}

/**
 * \brief Computes limited antidiffusive correction flux vector.
 *
 * This function computes the limited antidiffusive flux vector
 * \f$\bar{\mathbf{p}}\f$, which has the following entries:
 * \f[
 *   \bar{p}_i = \sum\limits_j L_{i,j} P_{i,j} \,.
 * \f]
 * The result is stored in \link flux_correction_vector \endlink.
 *
 * \pre This function assumes
 * - the antidiffusive correction flux matrix \f$\mathbf{P}\f$ has been computed
 *   and stored in \link flux_correction_matrix \endlink
 * - the limiter flux matrix \f$\mathbf{L}\f$ has been computed
 *   and stored in \link limiter_matrix \endlink
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
 * \brief Computes antidiffusive correction flux vector without any limitation.
 *
 * This function computes the full antidiffusive flux vector
 * \f$\mathbf{p}\f$, which has the following entries:
 * \f[
 *   p_i = \sum\limits_j P_{i,j} \,.
 * \f]
 * The result is stored in \link flux_correction_vector \endlink.
 *
 * \pre This function assumes that the antidiffusive correction flux matrix
 * \f$\mathbf{P}\f$ has been computed and stored in \link flux_correction_matrix
 * \endlink.
 */
template <int dim>
void FCT<dim>::compute_full_flux_correction_vector()
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
 * \brief Check to see if the FCT bounds were satisfied at all time steps.
 *
 * \param[in] new_solution  new solution vector \f$\mathbf{U}^{n+1}\f$
 *
 * \return flag telling whether imposed bounds were satisfied at all steps
 *
 * \note This function should not be used if characteristic limiting was
 *       used, since the imposed bounds are for the characteristic variables.
 */
template <int dim>
bool FCT<dim>::check_fct_bounds_satisfied(
  const Vector<double> & new_solution) const
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
      double value_i = new_solution(i);

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

  // return boolean for satisfaction of FCT bounds
  return fct_bounds_satisfied;
}

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
 * \brief Computes a local characteristic transformation matrix.
 *
 * This function throws an exception in this base class.
 *
 * \param[in] solution  solution vector \f$\mathbf{u}\f$ at which to evaluate
 *                      transformation
 *
 * \return transformation matrix \f$\mathbf{T}(\mathbf{u})\f$
 */
template <int dim>
FullMatrix<double> FCT<dim>::compute_transformation_matrix(
  const Vector<double> &) const
{
  // throw exception if not overridden by derived class
  AssertThrow(false, ExcNotImplemented());

  // return dummy matrix to avoid compiler warning about not returning anything
  FullMatrix<double> dummy_matrix(n_components, n_components);
  return dummy_matrix;
}

/**
 * \brief Computes a local characteristic transformation matrix inverse.
 *
 * This function throws an exception in this base class.
 *
 * \param[in] solution  solution vector \f$\mathbf{u}\f$ at which to evaluate
 *                      transformation
 *
 * \return transformation matrix inverse \f$\mathbf{T}^{-1}(\mathbf{u})\f$
 */
template <int dim>
FullMatrix<double> FCT<dim>::compute_transformation_matrix_inverse(
  const Vector<double> &) const
{
  // throw exception if not overridden by derived class
  AssertThrow(false, ExcNotImplemented());

  // return dummy matrix to avoid compiler warning about not returning anything
  FullMatrix<double> dummy_matrix(n_components, n_components);
  return dummy_matrix;
}

/**
 * \brief Transforms a vector in terms of conservative variables to a vector
 *        in terms of characteristic variables.
 *
 * This function applies a local transformation evaluated with the solution
 * \f$\mathbf{U}\f$ to characteristic variables on a vector \f$\mathbf{y}\f$:
 * \f[
 *   \hat{\mathbf{y}}_i = \mathbf{T}^{-1}(\mathbf{U}_i)\mathbf{y}_i \,,
 * \f]
 * where \f$\hat{\mathbf{y}}\f$ is the transformed vector, and \f$i\f$ is a
 * node index.
 *
 * \param[in] solution  solution vector \f$\mathbf{U}\f$ at which to evaluate
 *                      transformations
 * \param[in] vector_original  vector \f$\mathbf{y}\f$ to be transformed
 * \param[in] vector_transformed  transformed vector \f$\hat{\mathbf{y}}\f$
 */
template <int dim>
void FCT<dim>::transform_vector(const Vector<double> & solution,
                                const Vector<double> & vector_original,
                                Vector<double> & vector_transformed) const
{
  // loop over nodes
  for (unsigned int n = 0; n < n_nodes; ++n)
  {
    // assemble solution vector on node
    Vector<double> solution_node(n_components);
    for (unsigned int m = 0; m < n_components; ++m)
      solution_node[m] = solution[node_dof_indices[n][m]];

    // compute local transformation matrix inverse on node
    const FullMatrix<double> transformation_matrix_inverse =
      compute_transformation_matrix_inverse(solution_node);

    // apply local transformation matrix inverse to vector
    // loop over rows of matrix
    for (unsigned int m = 0; m < n_components; ++m)
    {
      // initialize sum for matrix-vector product for row
      const unsigned int i = node_dof_indices[n][m];
      vector_transformed[i] = 0.0;

      // loop over columns of matrix to compute matrix-vector product for row
      for (unsigned int mm = 0; mm < n_components; ++mm)
      {
        const unsigned int j = node_dof_indices[n][mm];
        vector_transformed[i] +=
          transformation_matrix_inverse[m][mm] * vector_original[j];
      }
    }
  }
}

/**
 * \brief Transforms a matrix in terms of conservative variables to a matrix
 *        in terms of characteristic variables.
 *
 * This function applies a local transformation evaluated with the solution
 * \f$\mathbf{U}\f$ to characteristic variables on a matrix \f$\mathbf{A}\f$:
 * \f[
 *   \hat{\mathbf{A}}_{i,j} = \mathbf{T}^{-1}(\mathbf{U}_i)\mathbf{A}_{i,j} \,,
 * \f]
 * where \f$\hat{\mathbf{A}}\f$ is the transformed matrix, and \f$i\f$ and
 * \f$j\f$ are node indices.
 *
 * \param[in] solution  solution vector \f$\mathbf{U}\f$ at which to evaluate
 *                      transformations
 * \param[in] matrix_original  matrix \f$\mathbf{A}\f$ to be transformed
 * \param[in] matrix_transformed  transformed matrix \f$\hat{\mathbf{A}}\f$
 */
template <int dim>
void FCT<dim>::transform_matrix(const Vector<double> & solution,
                                const SparseMatrix<double> & matrix_original,
                                SparseMatrix<double> & matrix_transformed) const
{
  // loop over nodes
  for (unsigned int n_i = 0; n_i < n_nodes; ++n_i)
  {
    // assemble solution vector on node
    Vector<double> solution_node(n_components);
    for (unsigned int m = 0; m < n_components; ++m)
      solution_node[m] = solution[node_dof_indices[n_i][m]];

    // compute local transformation matrix inverse on node
    const FullMatrix<double> transformation_matrix_inverse =
      compute_transformation_matrix_inverse(solution_node);

    // apply local transformation matrix inverse to vector
    // loop over components
    for (unsigned int m = 0; m < n_components; ++m)
    {
      // row index
      const unsigned int i = node_dof_indices[n_i][m];

      // loop over nonzero columns in row
      SparseMatrix<double>::const_iterator it_in = matrix_original.begin(i);
      SparseMatrix<double>::const_iterator it_end = matrix_original.end(i);
      SparseMatrix<double>::iterator it_out = matrix_transformed.begin(i);
      for (; it_in != it_end; ++it_in, ++it_out)
      {
        // get column index
        const unsigned int j = it_in->column();

        // get original matrix value for column
        const double value = it_in->value();

        // assert that if the value is greater than zero, then i and j correspond
        // to the same solution component
        Assert(std::abs(value) < 1.0e-15 || component_indices[j] == m,
               ExcInvalidState());

        // get node index corresponding to column
        const unsigned int n_j = node_indices[j];

        // compute transformed matrix value
        double transformed_value = 0.0;
        for (unsigned int l = 0; l < n_components; ++l)
        {
          const unsigned int k = node_dof_indices[n_i][l];
          const unsigned int p = node_dof_indices[n_j][l];
          const double value_original = matrix_original(k, p);

          transformed_value +=
            transformation_matrix_inverse[m][l] * value_original;
        }

        // store transformed matrix value
        it_out->value() = transformed_value;
      }
    }
  }
}

/**
 * \brief Creates lists of degree of freedom indices lists such as node list,
 *        component list, and nodal degrees of freedom.
 *
 * \c deal.II documentation states that \e local degree of freedom indices have a
 * a standard ordering: vertex 0 DoFs, vertex 1 DoFs, etc. Within each vertex,
 * it is assumed here that local ordering is by component index.
 */
template <int dim>
void FCT<dim>::create_dof_indices_lists()
{
  // initialize
  node_indices.resize(n_dofs);
  component_indices.resize(n_dofs);
  node_dof_indices.clear();

  // create vector for flagging whether a node index has been assigned or not
  std::vector<bool> assigned_node(n_dofs, false);

  // global degree of freedom indices in cell
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  unsigned int node_index = 0; // current node index
  Cell cell = dof_handler->begin_active();
  Cell endc = dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // get list of degrees of freedom on cell
    cell->get_dof_indices(local_dof_indices);

    // loop over nodes on cell
    for (unsigned int i_node = 0; i_node < dofs_per_cell_per_component; ++i_node)
    {
      // get DoF index of first DoF on node
      const unsigned int i_local0 = i_node * n_components;
      const unsigned int i0 = local_dof_indices[i_local0];

      // check if this DoF already has been assigned a node index
      if (!assigned_node[i0])
      {
        // vector of DoF indices on this node
        std::vector<unsigned int> this_node_dof_indices(n_components);

        // loop over DoFs on node
        for (unsigned int m = 0; m < n_components; ++m)
        {
          // determine global DoF index for component
          const unsigned int i_local = i_node * n_components + m;
          const unsigned int i = local_dof_indices[i_local];

          // assign node index for DoF
          node_indices[i] = node_index;

          // assign component index for DoF
          component_indices[i] = m;

          // assign DoF index to list of DoF indices on node
          this_node_dof_indices[m] = i;
        }

        // push back nodal DoF indices list
        node_dof_indices.push_back(this_node_dof_indices);

        // flag that node index has been assigned for the first DoF on node
        assigned_node[i0] = true;

        // increment node index
        node_index++;
      }
    }
  }

  // keep number of nodes
  n_nodes = node_index;
}

/**
 * \brief Checks conservation by ensuring that the sum of the limited flux
 *        correction sum is zero.
 *
 * \param[in] flux_vector limited correction flux sum vector
 */
template <int dim>
void FCT<dim>::check_conservation(const Vector<double> & flux_vector)
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
  limiter_matrix.print_formatted(limiter_out, 10, true, 0, "0", 1);
  limiter_out.close();
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
  postprocessor.output_dof_transient(solution_min,
                                     time,
                                     *dof_handler,
                                     "lower_fct_bound",
                                     lower_bound_component_names,
                                     bounds_transient_file_index,
                                     times_and_lower_bound_filenames,
                                     false,
                                     false);

  // output upper bound
  postprocessor.output_dof_transient(solution_max,
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
