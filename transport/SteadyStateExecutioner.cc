/**
 * Constructor.
 */
template<int dim>
SteadyStateExecutioner<dim>::SteadyStateExecutioner(
  const TransportParameters<dim> & parameters,
  const Triangulation<dim> & triangulation,
  const Tensor<1, dim> & transport_direction,
  const FunctionParser<dim> & cross_section_function,
  FunctionParser<dim> & source_function, Function<dim> & incoming_function,
  PostProcessor<dim> & postprocessor_) :
    Executioner<dim>(parameters, triangulation, transport_direction,
      cross_section_function, source_function, incoming_function, postprocessor_)
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
  Assert(this->parameters.viscosity_option != 2, ExcNotImplemented());
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

      // distribute constraints
      this->constraints.distribute(this->new_solution);

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

      // distribute constraints
      this->constraints.distribute(this->new_solution);

      break;
    }
    default :
    {
      ExcNotImplemented();
    }
  }

  // evaluate errors for convergence study
  this->postprocessor->evaluate_error(this->new_solution, this->dof_handler, *this->triangulation);

  // output grid and solution and print convergence results
  this->postprocessor->output_results(this->new_solution, this->dof_handler, *this->triangulation);
}
