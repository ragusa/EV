/**
 * \file SSPRKTimeIntegrator.cc
 * \brief Provides the function definitions for the SSPRKTimeIntegrator class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
SSPRKTimeIntegrator<dim>::SSPRKTimeIntegrator(
  const typename RunParameters<dim>::TemporalDiscretization & ssprk_method,
  const unsigned int & system_size,
  const LinearSolver<dim> & linear_solver,
  const SparsityPattern & sparsity_pattern)
  : n(system_size), linear_solver(linear_solver)
{
  // determine number of stages for chosen method
  switch (ssprk_method)
  {
    case RunParameters<dim>::TemporalDiscretization::FE:
      n_stages = 1;
      break;
    case RunParameters<dim>::TemporalDiscretization::SSP2:
      n_stages = 2;
      break;
    case RunParameters<dim>::TemporalDiscretization::SSP3:
      n_stages = 3;
      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // allocate memory for constants and stage solutions
  alpha.resize(n_stages);
  beta.resize(n_stages);
  c.resize(n_stages);
  u_stage.resize(n_stages + 1);
  for (unsigned int i = 0; i < n_stages + 1; ++i)
    u_stage[i].reinit(n);

  // assign RK parameters a, b, and c
  switch (ssprk_method)
  {
    case RunParameters<dim>::TemporalDiscretization::FE:
      alpha[0] = 0.0;
      beta[0] = 1.0;
      c[0] = 0.0;
      break;
    case RunParameters<dim>::TemporalDiscretization::SSP2:
      alpha[0] = 0.0;
      beta[0] = 1.0;
      c[0] = 0.0;
      alpha[1] = 0.5;
      beta[1] = 0.5;
      c[1] = 1.0;
      break;
    case RunParameters<dim>::TemporalDiscretization::SSP3:
      alpha[0] = 0.0;
      beta[0] = 1.0;
      c[0] = 0.0;
      alpha[1] = 0.75;
      beta[1] = 0.25;
      c[1] = 1.0;
      alpha[2] = 1.0 / 3.0;
      beta[2] = 2.0 / 3.0;
      c[2] = 0.5;
      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // initialize current stage index to 0
  current_stage = 0;

  // initialize time and time step size to zero
  t_old = 0.0;
  dt = 0.0;

  // allocate memory for temporary vectors
  intermediate_solution.reinit(n);
  system_rhs.reinit(n);

  // initialize sparse matrix
  system_matrix.reinit(sparsity_pattern);
}

/**
 * \brief Resets the current stage index, old solution, and time step size
 *
 * \param[in] old_solution old solution \f$\mathrm{U}^n\f$
 * \param[in] time_step_size time step size \f$\Delta t\f$
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::initialize_time_step(
  const Vector<double> & old_solution, const double & time_step_size)
{
  // update old time
  t_old += dt;

  // reset current stage index
  current_stage = 0;

  // set the old solution
  u_stage[0] = old_solution;

  // set the time step size
  dt = time_step_size;
}

/**
 * \brief Retrieves the current stage solution.
 *
 * \param[in] i stage index to retrieve
 * \param[out] solution vector for storage of retrieved stage solution
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::get_stage_solution(const unsigned int & i,
                                                  Vector<double> & solution) const
{
  // return stage solution
  solution = u_stage[i];
}

/**
 * \brief Retrieves the new solution.
 *
 * \param[out] solution vector for storage of new solution
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::get_new_solution(
  Vector<double> & new_solution) const
{
  // check that the correct number of stages were taken
  Assert(current_stage == n_stages, ExcInvalidState());

  // return final stage solution, which is new solution
  new_solution = u_stage[n_stages];
}

/**
 * \brief Advances one stage of a Runge-Kutta method.
 *
 * Here the following linear system is solved for the intermediate solution
 * \f$\tilde{\mathrm{U}}^i\f$:
 * \f[
 *   \mathrm{M}\tilde{\mathrm{U}}^i = \mathrm{M}\hat{\mathrm{U}}^{i-1} +
 *     \Delta t\left(\mathrm{r}(\hat{t}^i,\hat{\mathrm{U}}^{i-1})\right) \,.
 * \f]
 * If there are no additional steps (such as FCT) to take on the intermediate
 * solution, then the user specifies with a flag that the stage solution
 * can be completed here.
 *
 * \param[in] mass_matrix mass matrix
 * \param[in] ss_flux steady-state flux vector
 * \param[in] ss_rhs steady-state source vector
 * \param[in] diffusion_matrix diffusion matrix
 * \param[in] call_complete_stage_solution flag to call function to complete
 *            stage solution
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::step(const SparseMatrix<double> & mass_matrix,
                                    const Vector<double> & ss_flux,
                                    const Vector<double> & ss_rhs,
                                    const SparseMatrix<double> & diffusion_matrix,
                                    const bool & call_complete_stage_solution)
{
  // form transient rhs: system_rhs = M*u_old + dt*(ss_rhs - ss_flux - D*u_old)
  system_rhs = 0;
  mass_matrix.vmult(system_rhs, u_stage[current_stage]);
  system_rhs.add(dt, ss_rhs);
  system_rhs.add(-dt, ss_flux);
  diffusion_matrix.vmult(intermediate_solution, u_stage[current_stage]);
  system_rhs.add(-dt, intermediate_solution);

  // solve the linear system M*u_new = system_rhs
  system_matrix.copy_from(mass_matrix);
  double t_stage = get_stage_time();
  linear_solver.solve(
    system_matrix, intermediate_solution, system_rhs, true, t_stage);

  // compute new stage solution if there are no additional steps to be taken on
  // intermediate solution before the final combination
  if (call_complete_stage_solution)
    complete_stage_solution();
}

/**
 * \brief Computes stage solution \f$\hat{\mathrm{U}}^i\f$.
 *
 * At this stage, the intermediate stage solution \f$\tilde{\mathrm{U}}^i\f$
 * should have already been computed by solving the following linear system:
 * \f[
 *   \mathrm{M}\tilde{\mathrm{U}}^i = \mathrm{M}\hat{\mathrm{U}}^{i-1} +
 *     \Delta t\left(\mathrm{r}(\hat{t}^i,\hat{\mathrm{U}}^{i-1})\right) \,.
 * \f]
 * The stage solution is completed by performing the following linear
 * combination:
 * \f[
 *   \hat{\mathrm{U}}^i = \alpha_i\hat{\mathrm{U}}^0
 *     + \beta_i\tilde{\mathrm{U}^i} \,,
 * \f]
 * where \f$\hat{\mathrm{U}}^0\f$ is the old solution \f$\mathrm{U}^n\f$.
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::complete_stage_solution()
{
  // advance current stage index
  current_stage++;

  // perform combination
  u_stage[current_stage] = 0;
  u_stage[current_stage].add(alpha[current_stage - 1],
                             u_stage[0],
                             beta[current_stage - 1],
                             intermediate_solution);
}

/**
 * \brief Gets the intermediate stage solution \f$\tilde{\mathrm{U}}^i\f$.
 *
 * \param[out] solution vector for storage of intermediate stage solution.
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::get_intermediate_solution(
  Vector<double> & solution) const
{
  solution = intermediate_solution;
}

/**
 * \brief Sets the intermediate stage solution \f$\tilde{\mathrm{U}}^i\f$.
 *
 * \param[out] intermediate stage solution \f$\tilde{\mathrm{U}}^i\f$
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::set_intermediate_solution(
  const Vector<double> & solution)
{
  intermediate_solution = solution;
}

/**
 * \brief Gets the current stage time \f$\hat{t}^i = t^n + c_i\Delta t\f$.
 *
 * \return current stage time \f$\hat{t}^i = t^n + c_i\Delta t\f$
 */
template <int dim>
double SSPRKTimeIntegrator<dim>::get_stage_time() const
{
  return t_old + c[current_stage] * dt;
}

/**
 * \brief Sets a stage solution.
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::set_stage_solution(const unsigned int & i,
                                                  const Vector<double> & solution)
{
  u_stage[i] = solution;
}
