/** \brief constructor; initializes SSP Runge Kutta time integrator;
 *         assigns constant tables based on number of stages/order.
 */
template <int dim>
SSPRKTimeIntegrator<dim>::SSPRKTimeIntegrator(
  const SSPRKMethod & ssprk_method,
  const unsigned int & system_size,
  const LinearSolver<dim> & linear_solver,
  const SparsityPattern & sparsity_pattern)
  : n(system_size), linear_solver(linear_solver)
{
  // determine number of stages for chosen method
  switch (ssprk_method)
  {
    case SSPRKMethod::FE:
      n_stages = 1;
      break;
    case SSPRKMethod::SSP2:
      n_stages = 2;
      break;
    case SSPRKMethod::SSP3:
      n_stages = 3;
      break;
    default:
      Assert(false, ExcNotImplemented());
      break;
  }

  // allocate memory for constants and stage solutions
  a.reinit(n_stages);
  b.reinit(n_stages);
  c.reinit(n_stages);
  u_stage.resize(n_stages + 1);
  for (unsigned int i = 0; i < n_stages + 1; ++i)
    u_stage[i].reinit(n);

  // assign RK parameters a, b, and c
  switch (ssprk_method)
  {
    case SSPRKMethod::FE:
      a[0] = 0.0;
      b[0] = 1.0;
      c[0] = 0.0;
      break;
    case SSPRKMethod::SSP2:
      a[0] = 0.0;
      b[0] = 1.0;
      c[0] = 0.0;
      a[1] = 0.5;
      b[1] = 0.5;
      c[1] = 1.0;
      break;
    case SSPRKMethod::SSP3:
      a[0] = 0.0;
      b[0] = 1.0;
      c[0] = 0.0;
      a[1] = 0.75;
      b[1] = 0.25;
      c[1] = 1.0;
      a[2] = 1.0 / 3.0;
      b[2] = 2.0 / 3.0;
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

/** \brief resets the current stage index, old solution, and time step size
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

/** \brief Retrieves the current stage solution.
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::get_stage_solution(
  const unsigned int & i, Vector<double> & new_solution) const
{
  // return final stage solution, which is new solution
  new_solution = u_stage[i];
}

/** \brief retrieves the new solution.
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

/** \brief Advances one stage of a Runge-Kutta method.
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::step(const SparseMatrix<double> & mass_matrix,
                                    const SparseMatrix<double> & ss_matrix,
                                    const Vector<double> & ss_rhs,
                                    const bool & call_complete_stage_solution)
{
  // form transient rhs: system_rhs = M*u_old + dt*(ss_rhs - A*u_old)
  system_rhs = 0;
  system_rhs.add(dt, ss_rhs); //  now, system_rhs = dt*(ss_rhs)
  mass_matrix.vmult(intermediate_solution, u_stage[current_stage]);
  system_rhs.add(
    1.0, intermediate_solution); //  now, system_rhs = M*u_old + dt*(ss_rhs)
  ss_matrix.vmult(intermediate_solution, u_stage[current_stage]);
  system_rhs.add(
    -dt,
    intermediate_solution); //  now, system_rhs = M*u_old + dt*(ss_rhs - A*u_old)

  // solve the linear system M*u_new = system_rhs
  system_matrix.copy_from(mass_matrix);
  double t_stage = get_stage_time();
  linear_solver.solve(
    system_matrix, intermediate_solution, system_rhs, true, t_stage);

  // compute new stage solution if there are no additional steps to be taken on
  // intermediate solution
  if (call_complete_stage_solution)
    complete_stage_solution();
}

/** \brief Computes stage solution.
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::complete_stage_solution()
{
  // advance current stage index
  current_stage++;

  u_stage[current_stage] = 0;
  u_stage[current_stage].add(a[current_stage - 1],
                             u_stage[0],
                             b[current_stage - 1],
                             intermediate_solution);
}

/** \brief Gets the intermediate stage solution.
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::get_intermediate_solution(
  Vector<double> & solution) const
{
  solution = intermediate_solution;
}

/** \brief Sets the intermediate stage solution.
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::set_intermediate_solution(
  Vector<double> & solution)
{
  intermediate_solution = solution;
}

/** \brief Gets the current stage time.
 */
template <int dim>
double SSPRKTimeIntegrator<dim>::get_stage_time() const
{
  return t_old + c[current_stage] * dt;
}

/** \brief Sets a stage solution.
 */
template <int dim>
void SSPRKTimeIntegrator<dim>::set_stage_solution(const unsigned int & i,
                                                  const Vector<double> & solution)
{
  u_stage[i] = solution;
}