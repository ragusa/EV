/**
 * \brief Constructor.
 */
template <int dim>
TransportSSPRKExecutioner<dim>::TransportSSPRKExecutioner(
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
    ssprk(this->parameters.ssprk_discretization,
          this->n_dofs,
          this->linear_solver,
          this->constrained_sparsity_pattern)
{
  // create FCT object if needed
  if (this->parameters.scheme == Scheme::fct)
  {
    fct = std::make_shared<TransportExplicitEulerFCT<dim>>(
      this->parameters,
      *this->problem_parameters,
      this->dof_handler,
      this->fe,
      this->dirichlet_values,
      this->consistent_mass_matrix,
      this->lumped_mass_matrix);
  }

  // vector for old stage solution
  old_stage_solution.reinit(this->n_dofs);
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
void TransportSSPRKExecutioner<dim>::compute_new_solution(const double & dt,
                                                          const double & dt_old,
                                                          const double &,
                                                          const unsigned int & n)
{
  // initialize SSPRK time step
  ssprk.initialize_time_step(this->old_solution, dt);

  // loop over SSPRK stages to compute new solution
  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
  {
    // get stage time
    const double t_stage = ssprk.get_stage_time();

    // determine old stage solution and dt (needed for entropy viscosity)
    double old_stage_dt;
    if (i == 0)
    {
      old_stage_dt = dt_old;
    }
    else
    {
      ssprk.get_stage_solution(i - 1, this->older_solution);
      old_stage_dt = dt;
    }

    // get old stage solution
    ssprk.get_stage_solution(i, old_stage_solution);

    // compute steady-state rhs vector
    if (this->source_is_time_dependent)
      this->assembleSteadyStateRHS(this->ss_rhs, t_stage);

    // advance by an SSPRK step
    switch (this->parameters.scheme)
    {
      case Scheme::low:
        this->low_order_viscosity->update(
          old_stage_solution, this->older_solution, old_stage_dt, n);
        this->low_order_diffusion->compute_diffusion_matrix(
          old_stage_solution,
          this->low_order_viscosity,
          this->low_order_diffusion_matrix);
        ssprk.step(this->lumped_mass_matrix,
                   this->inviscid_ss_matrix,
                   this->low_order_diffusion_matrix,
                   this->ss_rhs,
                   true);
        break;
      case Scheme::high:
        this->high_order_viscosity->update(
          old_stage_solution, this->older_solution, old_stage_dt, n);
        this->high_order_diffusion->compute_diffusion_matrix(
          old_stage_solution,
          this->high_order_viscosity,
          this->high_order_diffusion_matrix);
        ssprk.step(this->consistent_mass_matrix,
                   this->inviscid_ss_matrix,
                   this->high_order_diffusion_matrix,
                   this->ss_rhs,
                   true);

        // add number of entropy viscosity iterations to total
        if (this->parameters.high_order_scheme == HighOrderScheme::entropy_visc)
          this->total_entropy_viscosity_iterations += ssprk.n_stages;

        break;
      case Scheme::fct:
        perform_fct_ssprk_step(dt, old_stage_dt, n);

        // add number of entropy viscosity iterations to total
        if (this->parameters.high_order_scheme == HighOrderScheme::entropy_visc)
          this->total_entropy_viscosity_iterations += ssprk.n_stages;

        break;
      default:
        AssertThrow(false, ExcNotImplemented());
        break;
    }
  }

  // retrieve the final solution
  ssprk.get_new_solution(this->new_solution);
}

/**
 * \brief Returns a pointer to FCT object.
 *
 * \return pointer to FCT object
 */
template <int dim>
std::shared_ptr<FCT<dim>> TransportSSPRKExecutioner<dim>::get_derived_fct() const
{
  return fct;
}

/**
 * \brief Performs an SSPRK step using FCT.
 *
 * \param[in] dt current time step size
 * \param[in] old_stage_dt time step size of previous SSPRK stage
 * \param[in] n time index
 */
template <int dim>
void TransportSSPRKExecutioner<dim>::perform_fct_ssprk_step(
  const double & dt, const double & old_stage_dt, const unsigned int & n)
{
  // update low-order diffusion
  this->low_order_viscosity->update(
    old_stage_solution, this->older_solution, old_stage_dt, n);
  this->low_order_diffusion->compute_diffusion_matrix(
    old_stage_solution,
    this->low_order_viscosity,
    this->low_order_diffusion_matrix);

  // update high-order diffusion
  this->high_order_viscosity->update(
    old_stage_solution, this->older_solution, old_stage_dt, n);
  this->high_order_diffusion->compute_diffusion_matrix(
    old_stage_solution,
    this->high_order_viscosity,
    this->high_order_diffusion_matrix);

  // compute high-order solution
  ssprk.step(this->consistent_mass_matrix,
             this->inviscid_ss_matrix,
             this->high_order_diffusion_matrix,
             this->ss_rhs,
             false);

  // get high-order solution
  ssprk.get_intermediate_solution(this->new_solution);

  // perform FCT
  // form rhs: system_rhs = M*u_old + dt*(ss_rhs - ss_flux - D*u_old + f)
  this->system_rhs = 0;
  this->lumped_mass_matrix.vmult(this->tmp_vector, this->old_solution);
  this->system_rhs.add(1.0, this->tmp_vector);
  this->system_rhs.add(dt, this->ss_rhs);
  this->inviscid_ss_matrix.vmult(this->inviscid_ss_flux, this->old_solution);
  this->system_rhs.add(-dt, this->inviscid_ss_flux);
  this->low_order_diffusion_matrix.vmult(this->tmp_vector, this->old_solution);
  this->system_rhs.add(-dt, this->tmp_vector);

  // compute antidiffusion vector
  Vector<double> & antidiffusion_vector = this->tmp_vector;
  fct->compute_antidiffusion_vector(this->new_solution,
                                    this->old_solution,
                                    dt,
                                    this->inviscid_ss_flux,
                                    this->reaction_vector,
                                    this->low_order_diffusion_matrix,
                                    this->high_order_diffusion_matrix,
                                    this->ss_rhs,
                                    antidiffusion_vector);

  // add antidiffusion vector
  this->system_rhs.add(dt, antidiffusion_vector);

  // solve the linear system M*u_new = system_rhs
  this->system_matrix.copy_from(this->lumped_mass_matrix);
  this->linear_solver.solve_with_dirichlet(
    this->system_matrix, this->new_solution, this->system_rhs);

  // set stage solution to be FCT solution for this stage
  ssprk.set_intermediate_solution(this->new_solution);

  // check FCT bounds
  if (this->parameters.check_fct_bounds)
    fct->check_bounds(this->new_solution);

  // finish computing stage solution
  ssprk.complete_stage_solution();
}
