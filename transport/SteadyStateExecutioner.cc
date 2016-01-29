/**
 * Constructor.
 */
template<int dim>
SteadyStateExecutioner<dim>::SteadyStateExecutioner(
  const TransportParameters<dim> & parameters,
  Triangulation<dim> & triangulation,
  const Tensor<1, dim> & transport_direction,
  const FunctionParser<dim> & cross_section_function,
  FunctionParser<dim> & source_function,
  Function<dim> & incoming_function,
  const double & domain_volume_,
  PostProcessor<dim> & postprocessor_) :
    Executioner<dim>(parameters, triangulation, transport_direction,
      cross_section_function, source_function, incoming_function, domain_volume_,
      postprocessor_)
{
}

/**
 * Destructor.
 */
template<int dim>
SteadyStateExecutioner<dim>::~SteadyStateExecutioner()
{
}

/**
 * Runs steady-state executioner.
 */
template<int dim>
void SteadyStateExecutioner<dim>::run()
{
  // entropy viscosity and FCT not yet implemented in steady-state
  Assert(this->parameters.viscosity_option != 3, ExcNotImplemented());
  Assert(this->parameters.viscosity_option != 4, ExcNotImplemented());

  // compute inviscid system matrix and steady-state right hand side (ss_rhs)
  this->assembleInviscidSteadyStateMatrix();
  this->assembleSteadyStateRHS(0.0);

  switch (this->parameters.viscosity_option)
  {
    case 0 :
    {
      // copy inviscid steady-state matrix to system matrix
      this->system_matrix.copy_from(this->inviscid_ss_matrix);

      // solve the linear system: ss_matrix*new_solution = ss_rhs
      this->linear_solver.solve(this->system_matrix, this->ss_rhs,
        this->new_solution);

      break;
    }
    case 1 :
    {
      // compute low-order viscosity
      LowOrderViscosity<dim> low_order_viscosity(this->n_cells,
        this->dofs_per_cell, this->dof_handler, this->constraints,
        this->inviscid_ss_matrix, this->low_order_diffusion_matrix,
        this->low_order_ss_matrix);

      // copy low-order steady-state matrix to system matrix
      this->system_matrix.copy_from(this->low_order_ss_matrix);

      // solve the linear system: ss_matrix*new_solution = ss_rhs
      this->linear_solver.solve(this->system_matrix, this->ss_rhs,
        this->new_solution);

      break;
    }
    case 2 :
    {
      // compute low-order viscosity
      LowOrderViscosity<dim> low_order_viscosity(
        this->n_cells, this->dofs_per_cell,
        this->dof_handler, this->constraints, this->inviscid_ss_matrix,
        this->low_order_diffusion_matrix, this->low_order_ss_matrix);

      // create entropy viscosity
      EntropyViscosity<dim> entropy_viscosity(
        this->fe, this->n_cells, this->dof_handler,
        this->constraints, this->cell_quadrature, this->face_quadrature,
        this->transport_direction, *this->cross_section_function,
        *this->source_function, this->parameters.entropy_string,
        this->parameters.entropy_derivative_string,
        this->parameters.entropy_residual_coefficient,
        this->parameters.jump_coefficient, this->domain_volume,
        this->parameters.EV_time_discretization, low_order_viscosity,
        this->inviscid_ss_matrix, this->high_order_diffusion_matrix,
        this->high_order_ss_matrix);

      // create nonlinear solver
      NonlinearSolver<dim> nonlinear_solver(
        this->parameters,
        this->constraints,
        this->dof_handler,
        *(this->incoming_function));

      // initialize guess for nonlinear solver
      this->new_solution = 0.0;
      nonlinear_solver.reset(this->new_solution);

      // begin iteration
      Vector<double> solution_change(this->n_dofs);
      Vector<double> tmp(this->n_dofs);
      bool converged = false;
      while (!converged)
      {
        // recompute high-order steady-state matrix A^(l)
        entropy_viscosity.recomputeHighOrderSteadyStateMatrix(
          this->new_solution);

        // create system matrix and rhs
        this->system_matrix.copy_from(this->high_order_ss_matrix);
        this->system_rhs = this->ss_rhs;

        // apply Dirichlet BC
        this->applyDirichletBC(
          this->system_matrix,
          this->system_rhs,
          this->new_solution);

        // compute linear residual: r^(l) = b - A^(l)*U^(l)
        this->system_matrix.vmult(tmp, this->new_solution);
        this->system_rhs.add(-1.0, tmp);

        // check convergence
        converged = nonlinear_solver.checkConvergence(this->system_rhs);

        // solve for solution update dU: A^(l)*dU = r^(l)
        this->linear_solver.solve(
          this->system_matrix,
          this->system_rhs,
          solution_change,
          false);

        // update solution: U^(l+1) = U^(l) + relax*dU
        this->new_solution.add(nonlinear_solver.relaxation_factor,
          solution_change);

        // distribute constraints
        this->constraints.distribute(this->new_solution);
      }

      break;
    }
    default :
    {
      ExcNotImplemented();
    }
  }

  // evaluate errors for convergence study
  this->postprocessor->evaluate_error(this->new_solution, this->dof_handler,
    *this->triangulation);

  // output grid and solution and print convergence results
  this->postprocessor->output_results(this->new_solution, this->dof_handler,
    *this->triangulation);
}
