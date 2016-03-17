/**
 * \brief Constructor.
 */
template <int dim>
TransportSteadyStateExecutioner<dim>::TransportSteadyStateExecutioner(
  const TransportRunParameters & parameters_,
  TransportProblemParameters<dim> & problem_parameters_,
  Triangulation<dim> & triangulation_,
  PostProcessor<dim> & postprocessor_)
  : TransportExecutioner<dim>(
      parameters_, problem_parameters_, triangulation_, postprocessor_)
{
}

/**
 * \brief Runs steady-state executioner.
 */
template <int dim>
void TransportSteadyStateExecutioner<dim>::run()
{
  // compute inviscid system matrix and steady-state right hand side (ss_rhs)
  this->assembleInviscidSteadyStateMatrix();
  this->assembleSteadyStateRHS(this->ss_rhs, 0.0);

  // compute low-order steady-state matrix since it is not solution-dependent
  // NOTE: arbitrary values are passed because the viscosity should not depend on
  // these parameters
  this->low_order_viscosity->update(
    this->new_solution, this->new_solution, 1.0, 1);
  this->low_order_diffusion->compute_diffusion_matrix(
    this->new_solution,
    this->low_order_viscosity,
    this->low_order_diffusion_matrix);
  this->low_order_ss_matrix.copy_from(this->inviscid_ss_matrix);
  this->low_order_ss_matrix.add(1.0, this->low_order_diffusion_matrix);

  switch (this->parameters.scheme)
  {
    case Scheme::low:
      compute_low_order_solution();
      break;
    case Scheme::high:
      compute_high_order_solution();
      break;
    case Scheme::fct:
      compute_high_order_solution();
      compute_fct_solution();
      break;
    default:
      throw ExcNotImplemented();
      break;
  }

  // evaluate errors for convergence study
  this->postprocessor->evaluate_error(
    this->new_solution, this->dof_handler, *this->triangulation);

  // output grid and solution and print convergence results
  this->postprocessor->output_results_if_last_cycle(
    this->new_solution, this->dof_handler, *this->triangulation);

  // print final solution if specified
  if (this->parameters.print_final_solution)
  {
    // set precision and format
    std::cout.precision(10);
    std::cout.setf(std::ios::scientific);

    // print each value of solution
    for (unsigned int j = 0; j < this->n_dofs; ++j)
      std::cout << this->new_solution[j] << std::endl;
  }
}

/**
 * \brief Computes the Galerkin solution.
 */
template <int dim>
void TransportSteadyStateExecutioner<dim>::compute_galerkin_solution()
{
  // copy inviscid steady-state matrix to system matrix
  this->system_matrix.copy_from(this->inviscid_ss_matrix);

  // copy steady-state right-hand-side vector to system rhs
  this->system_rhs = this->ss_rhs;

  // solve the linear system: ss_matrix*new_solution = ss_rhs
  this->linear_solver.solve_with_dirichlet(
    this->system_matrix, this->new_solution, this->system_rhs, true);
}

/**
 * \brief Computes the low-order solution.
 */
template <int dim>
void TransportSteadyStateExecutioner<dim>::compute_low_order_solution()
{
  // copy low-order steady-state matrix to system matrix
  this->system_matrix.copy_from(this->low_order_ss_matrix);

  // copy steady-state right-hand-side vector to system rhs
  this->system_rhs = this->ss_rhs;

  // solve the linear system: ss_matrix*new_solution = ss_rhs
  this->linear_solver.solve_with_dirichlet(
    this->system_matrix, this->new_solution, this->system_rhs, true);
}

/**
 * \brief Computes the high-order solution.
 */
template <int dim>
void TransportSteadyStateExecutioner<dim>::compute_high_order_solution()
{
  switch (this->parameters.high_order_scheme)
  {
    case HighOrderScheme::galerkin: // galerkin
      compute_galerkin_solution();
      break;
    case HighOrderScheme::entropy_visc: // entropy viscosity
      compute_entropy_viscosity_solution();
      break;
    default:
      throw ExcNotImplemented();
      break;
  }
}

/**
 * \brief Computes the entropy viscosity solution.
 */
template <int dim>
void TransportSteadyStateExecutioner<dim>::compute_entropy_viscosity_solution()
{
  // initialize guess for nonlinear solver
  this->new_solution = 0.0;
  this->nonlinear_solver.initialize(this->new_solution);

  // begin iteration
  bool converged = false;
  while (!converged)
  {
    // recompute high-order steady-state matrix A^(l)
    const unsigned int big_integer = 1e6; // dummy time step index
    this->high_order_viscosity->update(
      this->new_solution, this->new_solution, 1.0, big_integer);
    this->high_order_diffusion->compute_diffusion_matrix(
      this->new_solution,
      this->high_order_viscosity,
      this->high_order_diffusion_matrix);
    this->high_order_ss_matrix.copy_from(this->inviscid_ss_matrix);
    this->high_order_ss_matrix.add(1.0, this->high_order_diffusion_matrix);

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
void TransportSteadyStateExecutioner<dim>::compute_fct_solution()
{
  /*
    // create FCT object
    FCT<dim> fct(this->dof_handler,
                 *this->triangulation,
                 this->linear_solver,
                 this->constrained_sparsity_pattern,
                 this->dirichlet_nodes,
                 this->n_dofs,
                 this->dofs_per_cell,
                 false,
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
    if (this->parameters.output_fct_bounds)
      fct.output_bounds(*(this->postprocessor));

    // check FCT bounds
    fct.check_fct_bounds(this->new_solution);
  */
}
