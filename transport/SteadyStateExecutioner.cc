/**
 * Constructor.
 */
template <int dim>
SteadyStateExecutioner<dim>::SteadyStateExecutioner(
  const TransportParameters<dim> & parameters,
  Triangulation<dim> & triangulation,
  const Tensor<1, dim> & transport_direction,
  const FunctionParser<dim> & cross_section_function,
  FunctionParser<dim> & source_function,
  Function<dim> & incoming_function,
  const double & domain_volume_,
  PostProcessor<dim> & postprocessor_)
  : Executioner<dim>(parameters,
                     triangulation,
                     transport_direction,
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

  // solve the linear system: ss_matrix*new_solution = ss_rhs
  this->linear_solver.solve(
    this->system_matrix, this->new_solution, this->ss_rhs, true);
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
               this->parameters.do_not_limit);

  // compute flux corrections
  fct.compute_flux_corrections_ss(this->new_solution,
                                  this->low_order_diffusion_matrix,
                                  this->high_order_diffusion_matrix);

  // copy low-order steady-state matrix to system matrix
  this->system_matrix.copy_from(this->low_order_ss_matrix);

  // initialize guess for nonlinear solver
  this->new_solution = 0.0;
  this->nonlinear_solver.initialize(this->new_solution);

  // begin iteration
  bool converged = false;
  while (!converged)
  {
    // compute max principle min and max values
    fct.compute_bounds_ss(
      this->new_solution, this->low_order_ss_matrix, this->ss_rhs);

    // compute limited flux correction sum and add it to rhs
    fct.compute_limiting_coefficients_ss(
      this->new_solution, this->low_order_ss_matrix, this->ss_rhs);

    // create system rhs
    this->system_rhs = this->ss_rhs;
    this->system_rhs.add(1.0, fct.get_flux_correction_vector());

    // apply Dirichlet BC here
    this->applyDirichletBC(
      this->system_matrix, this->system_rhs, this->new_solution);

    // check convergence and perform update if necessary
    converged =
      this->nonlinear_solver.update(this->system_matrix, this->system_rhs);
  }
}
