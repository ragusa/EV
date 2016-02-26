/**
 * Constructor.
 */
template <int dim>
SteadyStateExecutioner<dim>::SteadyStateExecutioner(
  const TransportParameters<dim> & parameters,
  Triangulation<dim> & triangulation,
  const Tensor<1, dim> & transport_direction,
  const double & transport_speed,
  const FunctionParser<dim> & cross_section_function,
  FunctionParser<dim> & source_function,
  Function<dim> & incoming_function,
  const double & domain_volume_,
  PostProcessor<dim> & postprocessor_)
  : Executioner<dim>(parameters,
                     triangulation,
                     transport_direction,
                     transport_speed,
                     cross_section_function,
                     source_function,
                     incoming_function,
                     domain_volume_,
                     postprocessor_)
{
}

/**
 * Runs steady-state executioner.
 */
template <int dim>
void SteadyStateExecutioner<dim>::run()
{
  // compute inviscid system matrix and steady-state right hand side (ss_rhs)
  this->assembleInviscidSteadyStateMatrix();
  this->assembleSteadyStateRHS(this->ss_rhs, 0.0);

  switch (this->parameters.viscosity_option)
  {
    case 0: // Galerkin
    {
      // compute Galerkin solution
      compute_galerkin_solution();

      break;
    }
    case 1: // Low-Order
    {
      // compute low-order viscosity
      LowOrderViscosity<dim> low_order_viscosity(this->n_cells,
                                                 this->dofs_per_cell,
                                                 this->dof_handler,
                                                 this->constraints,
                                                 this->inviscid_ss_matrix,
                                                 this->low_order_diffusion_matrix,
                                                 this->low_order_ss_matrix);

      // compute the low-order solution
      compute_low_order_solution();

      break;
    }
    case 2: // Entropy Viscosity
    {
      // compute entropy viscosity solution
      compute_entropy_viscosity_solution();

      break;
    }
    case 3: // EV-FCT
    {
      // compute entropy viscosity solution
      compute_entropy_viscosity_solution();

      // compute FCT solution
      compute_FCT_solution();

      break;
    }
    case 4: // Galerkin FCT
    {
      // compute Galerkin solution
      compute_galerkin_solution();

      // compute low-order viscosity
      LowOrderViscosity<dim> low_order_viscosity(this->n_cells,
                                                 this->dofs_per_cell,
                                                 this->dof_handler,
                                                 this->constraints,
                                                 this->inviscid_ss_matrix,
                                                 this->low_order_diffusion_matrix,
                                                 this->low_order_ss_matrix);

      // compute FCT solution
      compute_FCT_solution();

      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
    }
  }

  // evaluate errors for convergence study
  this->postprocessor->evaluate_error(
    this->new_solution, this->dof_handler, *this->triangulation);

  // output grid and solution and print convergence results
  this->postprocessor->output_results(
    this->new_solution, this->dof_handler, *this->triangulation);
}

/**
 * \brief Computes the Galerkin solution.
 */
template <int dim>
void SteadyStateExecutioner<dim>::compute_galerkin_solution()
{
  // copy inviscid steady-state matrix to system matrix
  this->system_matrix.copy_from(this->inviscid_ss_matrix);

  // copy steady-state right-hand-side vector to system rhs
  this->system_rhs = this->ss_rhs;

  // solve the linear system: ss_matrix*new_solution = ss_rhs
  this->linear_solver.solve(
    this->system_matrix, this->new_solution, this->system_rhs, true);
}

/**
 * \brief Computes the low-order solution.
 */
template <int dim>
void SteadyStateExecutioner<dim>::compute_low_order_solution()
{
  // copy low-order steady-state matrix to system matrix
  this->system_matrix.copy_from(this->low_order_ss_matrix);

  // solve the linear system: ss_matrix*new_solution = ss_rhs
  this->linear_solver.solve(
    this->system_matrix, this->new_solution, this->ss_rhs, true);
}

/**
 * \brief Computes the entropy viscosity solution.
 */
template <int dim>
void SteadyStateExecutioner<dim>::compute_entropy_viscosity_solution()
{
  // compute low-order viscosity
  LowOrderViscosity<dim> low_order_viscosity(this->n_cells,
                                             this->dofs_per_cell,
                                             this->dof_handler,
                                             this->constraints,
                                             this->inviscid_ss_matrix,
                                             this->low_order_diffusion_matrix,
                                             this->low_order_ss_matrix);

  // create entropy viscosity
  EntropyViscosity<dim> entropy_viscosity(
    this->fe,
    this->n_cells,
    this->dof_handler,
    this->constraints,
    this->cell_quadrature,
    this->face_quadrature,
    this->transport_direction,
    this->transport_speed,
    *this->cross_section_function,
    *this->source_function,
    this->parameters.entropy_string,
    this->parameters.entropy_derivative_string,
    this->parameters.entropy_residual_coefficient,
    this->parameters.jump_coefficient,
    this->domain_volume,
    this->parameters.entropy_temporal_discretization,
    low_order_viscosity,
    this->inviscid_ss_matrix,
    this->high_order_diffusion_matrix,
    this->high_order_ss_matrix);

  // initialize guess for nonlinear solver
  this->new_solution = 0.0;
  this->nonlinear_solver.initialize(this->new_solution);

  // begin iteration
  bool converged = false;
  while (!converged)
  {
    // recompute high-order steady-state matrix A^(l)
    entropy_viscosity.recompute_high_order_ss_matrix(this->new_solution);

    // create system matrix and rhs
    this->system_matrix.copy_from(this->high_order_ss_matrix);
    this->system_rhs = this->ss_rhs;

    // apply Dirichlet BC
    this->applyDirichletBC(
      this->system_matrix, this->system_rhs, this->new_solution);

    // check convergence and perform update if necessary
    converged =
      this->nonlinear_solver.update(this->system_matrix, this->system_rhs);
  }
}

/**
 * \brief Solves the steady-state FCT system.
 */
template <int dim>
void SteadyStateExecutioner<dim>::compute_FCT_solution()
{
  // create FCT object
  FCT<dim> fct(this->dof_handler,
               *this->triangulation,
               this->linear_solver,
               this->constrained_sparsity_pattern,
               this->dirichlet_nodes,
               this->n_dofs,
               this->dofs_per_cell,
               this->parameters.do_not_limit,
               this->parameters.include_analytic_bounds,
               this->fe,
               this->cell_quadrature,
               *this->cross_section_function,
               *this->source_function);

  // check if high-order solution satisfies bounds - if so, do not use FCT
  bool skip_fct = false;
  if (this->parameters.skip_fct_if_bounds_satisfied)
  {
    // compute max principle min and max values
    fct.compute_bounds_ss(
      this->new_solution, this->low_order_ss_matrix, this->ss_rhs);

    // check bounds
    skip_fct = fct.check_fct_bounds(this->new_solution);
  }

  if (!skip_fct)
  {
    // compute flux corrections
    fct.compute_flux_corrections_ss(this->new_solution,
                                    this->low_order_diffusion_matrix,
                                    this->high_order_diffusion_matrix);

    // initialize guess for nonlinear solver
    switch (this->parameters.fct_initialization_option)
    {
      case FCTInitializationOption::zero:
      {
        // set to zero
        this->new_solution = 0.0;

        break;
      }
      case FCTInitializationOption::low:
      {
        // compute low-order solution
        compute_low_order_solution();
        break;
      }
      case FCTInitializationOption::high:
      {
        // do nothing, solution vector already contains high-order solution
        break;
      }
      default:
      {
        Assert(false, ExcNotImplemented());
      }
    }
    this->nonlinear_solver.initialize(this->new_solution);

    // initialize cumulative antidiffusion vector
    this->cumulative_antidiffusion = 0.0;

    // begin iteration
    bool converged = false;
    while (!converged)
    {
      // compute max principle min and max values: W^(l)
      fct.compute_bounds_ss(
        this->new_solution, this->low_order_ss_matrix, this->ss_rhs);

      // compute limited flux bounds: Q^(l)
      fct.compute_limited_flux_bounds_ss(
        this->new_solution, this->low_order_ss_matrix, this->ss_rhs);

      // compute limited flux correction matrix and sum: dp^(l) and dP^(l)
      fct.compute_limited_fluxes();

      // if using cumulative antidiffusion algorithm, then update cumulative
      // antidiffusion and remainder antidiffusive fluxes
      if (this->parameters.use_cumulative_antidiffusion_algorithm)
      {
        // add to cumulative antidiffusion vector: p^(l+1) = p^(l) + dp^(l)
        this->cumulative_antidiffusion.add(1.0, fct.get_limited_flux_vector());

        // subtract used antidiffusive flux: DP^(l+1) = DP^(l) - dP^(l)
        fct.subtract_limited_flux_correction_matrix();
      }
      else
      {
        // throw away accumulated antidiffusive flux
        this->cumulative_antidiffusion = fct.get_limited_flux_vector();
      }

      // create system rhs: s^(l) = b + p^(l+1)
      this->system_rhs = this->ss_rhs;
      this->system_rhs.add(1.0, this->cumulative_antidiffusion);

      // create system matrix. Note that although the underlying system matrix
      // does not change in each iteration, the Dirichlet-modified system matrix
      // does change
      this->system_matrix.copy_from(this->low_order_ss_matrix);

      // apply Dirichlet BC here
      this->applyDirichletBC(
        this->system_matrix, this->system_rhs, this->new_solution);

      // check convergence and perform update if necessary
      converged =
        this->nonlinear_solver.update(this->system_matrix, this->system_rhs);
    }
  }

  // output FCT bounds if requested
  if (this->parameters.output_DMP_bounds)
    fct.output_bounds(*(this->postprocessor));

  // check FCT bounds
  fct.check_fct_bounds(this->new_solution);
}
