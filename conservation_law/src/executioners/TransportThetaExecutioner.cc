/**
 * \brief Constructor.
 */
template <int dim>
TransportThetaExecutioner<dim>::TransportThetaExecutioner(
  const TransportRunParameters & parameters_,
  TransportProblemParameters<dim> & problem_parameters_,
  Triangulation<dim> & triangulation_,
  PostProcessor<dim> & postprocessor_,
  const double & nominal_dt_)
  : TransportTransientExecutioner<dim>(parameters_,
                                       problem_parameters_,
                                       triangulation_,
                                       postprocessor_,
                                       nominal_dt_),
    theta(parameters_.theta)
{
  // create FCT object if needed
  if (this->parameters.scheme == Scheme::fct)
  {
    fct = std::make_shared<TransportThetaFCT<dim>>(this->parameters,
                                                   *this->problem_parameters,
                                                   this->dof_handler,
                                                   this->fe,
                                                   this->dirichlet_values,
                                                   this->consistent_mass_matrix,
                                                   this->lumped_mass_matrix);
  }

  // initialize matrices
  high_order_diffusion_matrix_new.reinit(this->constrained_sparsity_pattern);
  high_order_ss_matrix_new.reinit(this->constrained_sparsity_pattern);
}

/**
 * \brief Computes new solution for a time step.
 *
 * \param[in] dt  current time step size \f$\Delta t^n\f$
 * \param[in] dt_old  old time step size \f$\Delta t^{n-1}\f$
 * \param[in] t_old  old time \f$t^n\f$
 * \param[in] n  time index
 */
template <int dim>
void TransportThetaExecutioner<dim>::compute_new_solution(const double & dt,
                                                          const double &,
                                                          const double & t_old,
                                                          const unsigned int & n)
{
  // compute new time
  const double t_new = t_old + dt;

  // compute new solution
  switch (this->parameters.scheme)
  {
    case Scheme::low:
      compute_low_order_solution(dt, t_new);
      break;
    case Scheme::high:
      compute_high_order_solution(dt, t_new, n);
      break;
    case Scheme::fct:
      compute_high_order_solution(dt, t_new, n);
      compute_fct_solution(dt, t_old);
      break;
    default:
      AssertThrow(false, ExcNotImplemented());
      break;
  }

  // save old quantities
  this->ss_rhs = this->ss_rhs_new;
  this->high_order_diffusion_matrix.copy_from(high_order_diffusion_matrix_new);
  this->high_order_ss_matrix.copy_from(high_order_ss_matrix_new);
}

/**
 * \brief Returns a pointer to FCT object.
 *
 * \return pointer to FCT object
 */
template <int dim>
std::shared_ptr<FCT<dim>> TransportThetaExecutioner<dim>::get_derived_fct() const
{
  return fct;
}

/**
 * \brief Takes time step with the Galerkin method using a theta method.
 *
 * \param[in] dt time step size
 * \param[in] t_new new time
 */
template <int dim>
void TransportThetaExecutioner<dim>::compute_galerkin_solution(
  const double & dt, const double & t_new)
{
  // compute new steady-state right-hand-side vector b_new
  if (this->source_is_time_dependent)
    this->assembleSteadyStateRHS(this->ss_rhs_new, t_new);

  // compute system matrix
  // matrix = M^C + theta*dt*A
  this->system_matrix.copy_from(this->consistent_mass_matrix);
  this->system_matrix.add(theta * dt, this->inviscid_ss_matrix);

  // compute system right-hand-side vector
  // rhs = M^C*u_old + dt*((1-theta)*b_old + theta*b_new - (1-theta)*A*u_old)
  this->consistent_mass_matrix.vmult(this->system_rhs, this->old_solution);
  this->inviscid_ss_matrix.vmult(this->tmp_vector, this->old_solution);
  this->system_rhs.add(-(1.0 - theta) * dt, this->tmp_vector);
  this->system_rhs.add((1.0 - theta) * dt, this->ss_rhs);
  this->system_rhs.add(theta * dt, this->ss_rhs_new);

  // solve linear system
  this->linear_solver.solve_with_dirichlet(
    this->system_matrix, this->new_solution, this->system_rhs, true, t_new);
}

/**
 * \brief Takes time step with the low-order method using a theta method.
 *
 * \param[in] dt     new time step size \f$\Delta t^n\equiv t^{n+1}-t^n\f$
 * \param[in] t_new  new time \f$t^{n+1}\f$
 */
template <int dim>
void TransportThetaExecutioner<dim>::compute_low_order_solution(
  const double & dt, const double & t_new)
{
  // compute new steady-state right-hand-side vector b_new
  if (this->source_is_time_dependent)
    this->assembleSteadyStateRHS(this->ss_rhs_new, t_new);

  // compute system matrix
  // matrix = M^L + theta*dt*A^L
  this->system_matrix.copy_from(this->lumped_mass_matrix);
  this->system_matrix.add(theta * dt, this->low_order_ss_matrix);

  // compute system right-hand-side vector
  // rhs = M^L*u_old + dt*((1-theta)*b_old + theta*b_new - (1-theta)*A^L*u_old)
  this->lumped_mass_matrix.vmult(this->system_rhs, this->old_solution);
  this->low_order_ss_matrix.vmult(this->tmp_vector, this->old_solution);
  this->system_rhs.add(-(1.0 - theta) * dt, this->tmp_vector);
  this->system_rhs.add((1.0 - theta) * dt, this->ss_rhs);
  this->system_rhs.add(theta * dt, this->ss_rhs_new);

  // solve linear system
  this->linear_solver.solve_with_dirichlet(
    this->system_matrix, this->new_solution, this->system_rhs, true, t_new);
}

/**
 * \brief Takes time step with a high-order method using a theta method.
 *
 * \param[in] dt     new time step size \f$\Delta t^n\equiv t^{n+1}-t^n\f$
 * \param[in] t_new  new time \f$t^{n+1}\f$
 * \param[in] n      time step index
 */
template <int dim>
void TransportThetaExecutioner<dim>::compute_high_order_solution(
  const double & dt, const double & t_new, const unsigned int & n)
{
  switch (this->parameters.high_order_scheme)
  {
    case HighOrderScheme::galerkin: // galerkin
      compute_galerkin_solution(dt, t_new);
      break;
    case HighOrderScheme::entropy_visc: // entropy viscosity
      compute_entropy_viscosity_solution(dt, t_new, n);
      break;
    default:
      AssertThrow(false, ExcNotImplemented());
      break;
  }
}

/**
 * \brief Takes time step with the entropy viscosity method using a theta method.
 *
 * \param[in] dt     new time step size \f$\Delta t^n    \equiv t^{n+1}-t^n\f$
 * \param[in] t_new  new time \f$t^{n+1}\f$
 * \param[in] n      time step index
 */
template <int dim>
void TransportThetaExecutioner<dim>::compute_entropy_viscosity_solution(
  const double & dt, const double & t_new, const unsigned int & n)
{
  // references
  Vector<double> & ss_rhs_old = this->ss_rhs;
  SparseMatrix<double> & high_order_ss_matrix_old = this->high_order_ss_matrix;

  // compute new steady-state right-hand-side vector b_new
  if (this->source_is_time_dependent)
    this->assembleSteadyStateRHS(this->ss_rhs_new, t_new);

  // initialize guess for nonlinear solver
  this->new_solution = 0.0;
  this->nonlinear_solver.initialize(this->new_solution);

  // begin iteration
  bool converged = false;
  while (!converged)
  {
    // compute new entropy viscosity and high order diffusion matrix
    this->high_order_viscosity->update(
      this->new_solution, this->old_solution, dt, n + 1);
    this->high_order_diffusion->compute_diffusion_matrix(
      this->new_solution,
      this->high_order_viscosity,
      high_order_diffusion_matrix_new);

    // compute new high-order steady-state matrix
    high_order_ss_matrix_new.copy_from(this->inviscid_ss_matrix);
    high_order_ss_matrix_new.add(1.0, high_order_diffusion_matrix_new);

    // compute system matrix
    // matrix = M^C + theta*dt*A^H
    this->system_matrix.copy_from(this->consistent_mass_matrix);
    this->system_matrix.add(theta * dt, high_order_ss_matrix_new);

    // compute system right-hand-side vector
    // rhs = M^C*u_old + dt*((1-theta)*b_old + theta*b_new - (1-theta)*A^H*u_old)
    this->consistent_mass_matrix.vmult(this->system_rhs, this->old_solution);
    high_order_ss_matrix_old.vmult(this->tmp_vector, this->old_solution);
    this->system_rhs.add(-(1.0 - theta) * dt, this->tmp_vector);
    this->system_rhs.add((1.0 - theta) * dt, ss_rhs_old);
    this->system_rhs.add(theta * dt, this->ss_rhs_new);

    // apply Dirichlet BC
    this->applyDirichletBC(
      this->system_matrix, this->system_rhs, this->new_solution);

    // check convergence and perform update if necessary
    converged =
      this->nonlinear_solver.update(this->system_matrix, this->system_rhs);
  }
}

/**
 * \brief Computes the FCT solution using a theta method.
 *
 * \param[in] dt  time step size
 */
template <int dim>
void TransportThetaExecutioner<dim>::compute_fct_solution(const double & dt,
                                                          const double &)
{
  // references
  const Vector<double> & ss_rhs_old = this->ss_rhs;
  const SparseMatrix<double> & high_order_diffusion_matrix_old =
    this->high_order_diffusion_matrix;

  // compute flux corrections
  fct->compute_antidiffusion_matrix(this->new_solution,
                                    this->old_solution,
                                    dt,
                                    this->low_order_diffusion_matrix,
                                    this->low_order_diffusion_matrix,
                                    high_order_diffusion_matrix_old,
                                    high_order_diffusion_matrix_new);

  // compute system matrix
  this->system_matrix.copy_from(this->lumped_mass_matrix);
  this->system_matrix.add(theta * dt, this->low_order_ss_matrix);

  // initialize guess for nonlinear solver
  this->new_solution = 0.0;
  this->nonlinear_solver.initialize(this->new_solution);

  // begin iteration
  bool converged = false;
  while (!converged)
  {
    // compute antidiffusion vector
    fct->compute_antidiffusion_vector(this->new_solution,
                                      this->old_solution,
                                      dt,
                                      this->low_order_ss_matrix,
                                      this->ss_rhs_new,
                                      ss_rhs_old,
                                      this->antidiffusion_vector);

    // create system rhs
    this->system_rhs = 0;
    this->system_rhs.add((1.0 - theta) * dt, ss_rhs_old);
    this->system_rhs.add(theta * dt, this->ss_rhs_new);
    this->lumped_mass_matrix.vmult(this->tmp_vector, this->old_solution);
    this->system_rhs.add(1.0, this->tmp_vector);
    this->low_order_ss_matrix.vmult(this->tmp_vector, this->old_solution);
    this->system_rhs.add(-(1.0 - theta) * dt, this->tmp_vector);
    this->system_rhs.add(dt, this->antidiffusion_vector);

    // create system matrix. Note that although the underlying system matrix
    // does not change in each iteration, the Dirichlet-modified system matrix
    // does change
    this->system_matrix.copy_from(this->lumped_mass_matrix);
    this->system_matrix.add(theta * dt, this->low_order_ss_matrix);

    // apply Dirichlet BC here
    this->applyDirichletBC(
      this->system_matrix, this->system_rhs, this->new_solution);

    // check convergence and perform update if necessary
    converged =
      this->nonlinear_solver.update(this->system_matrix, this->system_rhs);
  }

  // check FCT solution
  if (this->parameters.check_fct_bounds)
    fct->check_bounds(this->new_solution);
}
