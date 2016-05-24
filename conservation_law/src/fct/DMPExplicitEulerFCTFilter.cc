/**
 * \file DMPExplicitEulerFCTFilter.cc
 * \brief Provides the function definitions for the DMPExplicitEulerFCTFilter
 * class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] run_parameters_  run parameters
 * \param[in] limiter_  limiter
 * \param[in] dof_handler_  degree of freedom handler
 * \param[in] fe_  finite element system
 * \param[in] lumped_mass_matrix_  lumped mass matrix \f$\mathbf{M}^L\f$
 * \param[in] dirichlet_values_  map of DoF indices to Dirichlet values
 */
template <int dim>
DMPExplicitEulerFCTFilter<dim>::DMPExplicitEulerFCTFilter(
  const RunParameters & run_parameters_,
  const std::shared_ptr<Limiter<dim>> limiter_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const SparseMatrix<double> & lumped_mass_matrix_,
  const std::map<unsigned int, double> & dirichlet_values_)
  : ExplicitEulerFCTFilter<dim>(run_parameters_,
                                limiter_,
                                dof_handler_,
                                fe_,
                                lumped_mass_matrix_,
                                dirichlet_values_)
{
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
 * \param[in] dt  time step size \f$\Delta t\f$
 * \param[in] ss_reaction  steady-state reaction vector \f$\mathbf{\sigma}\f$,
 *            where \f$\sigma_i = \int\limits_{S_i}
 *            \varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *            = \sum\limits_j\int\limits_{S_{i,j}}
 *            \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \f$
 * \param[in] ss_rhs  old steady-state right hand side vector \f$\mathbf{b}^n\f$
 * \param[in] t_old  old time
 */
template <int dim>
void DMPExplicitEulerFCTFilter<dim>::compute_solution_bounds(
  const Vector<double> & old_solution,
  const double & dt,
  const Vector<double> & ss_reaction,
  const Vector<double> & ss_rhs,
  const double &)
{
  // create references
  Vector<double> & solution_min = this->solution_bounds.lower;
  Vector<double> & solution_max = this->solution_bounds.upper;
  const SparseMatrix<double> & lumped_mass_matrix = *this->lumped_mass_matrix;

  // compute minimum and maximum values of solution
  this->compute_min_and_max_of_dof_vector(
    old_solution, solution_min, solution_max);

  /*
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
  */

  // At this point, the min/max values of the old solution in the support
  // of test function i are stored in solution_min(i) and solution_max(i).
  // Now these values are multiplied by (1-dt/m(i))*sigma(i) and
  // added to dt/m(i)*b(i).

  // compute the upper and lower bounds for the FCT solution
  for (unsigned int i = 0; i < this->n_dofs; ++i)
  {
    // if not a Dirichlet value
    if (this->dirichlet_values->find(i) == this->dirichlet_values->end())
    {
      solution_max(i) =
        solution_max(i) * (1.0 - dt / lumped_mass_matrix(i, i) * ss_reaction(i)) +
        dt / lumped_mass_matrix(i, i) * ss_rhs(i);

      solution_min(i) =
        solution_min(i) * (1.0 - dt / lumped_mass_matrix(i, i) * ss_reaction(i)) +
        dt / lumped_mass_matrix(i, i) * ss_rhs(i);
    }
    else
    {
      // get Dirichlet value
      const double value = (*this->dirichlet_values).at(i);
      solution_max(i) = value;
      solution_min(i) = value;
    }
  }
}
