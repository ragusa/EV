template<int dim>
FCT<dim>::FCT(
   const SparseMatrix<double> &lumped_mass_matrix)
{
}

template<int dim>
FCT<dim>::~FCT()
{
}

template<int dim>
FCT<dim>::solve_FCT_system(const Vector<double> &old_solution,
                           const SparseMatrix<double> &low_order_ss_matrix,
                           const Vector<double> &ss_rhs,
                           const double         &dt)
{
   // form transient rhs: system_rhs = M*u_old + dt*(ss_rhs - (A+D)*u_old)
   system_rhs = 0;
   system_rhs.add(dt, ss_rhs); //       now, system_rhs = dt*(ss_rhs)
   lumped_mass_matrix.vmult(tmp_vector, old_solution);
   system_rhs.add(1.0, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs)
   low_order_ss_matrix.vmult(tmp_vector, old_solution);
   system_rhs.add(-dt, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs - A*u_old)

   // compute flux corrections
   compute_flux_corrections(dt);

   // compute max principle min and max values
   compute_bounds(dt);

   // compute Q+
   Q_plus = 0;
   lumped_mass_matrix.vmult(tmp_vector, solution_max);
   Q_plus.add(1.0/dt, tmp_vector);
   lumped_mass_matrix.vmult(tmp_vector, old_solution);
   Q_plus.add(-1.0/dt,tmp_vector);
   low_order_ss_matrix.vmult(tmp_vector, old_solution);
   Q_plus.add(1.0,tmp_vector);
   Q_plus.add(-1.0,ss_rhs);
   
   // compute Q-
   Q_minus = 0;
   lumped_mass_matrix.vmult(tmp_vector, solution_max);
   Q_minus.add(1.0/dt, tmp_vector);
   lumped_mass_matrix.vmult(tmp_vector, old_solution);
   Q_minus.add(-1.0/dt,tmp_vector);
   low_order_ss_matrix.vmult(tmp_vector, old_solution);
   Q_minus.add(1.0,tmp_vector);
   Q_minus.add(-1.0,ss_rhs);

   // compute limited flux correction sum and add it to rhs
   compute_limiting_coefficients();
   system_rhs.add(dt, flux_correction_vector);   // now, system_rhs = M*u_old + dt*(ss_rhs - A*u_old + f)

   // solve the linear system M*u_new = system_rhs
   system_matrix.copy_from(lumped_mass_matrix);
   apply_Dirichlet_BC(system_matrix, new_solution, system_rhs); // apply Dirichlet BC
   solve_linear_system(system_matrix, system_rhs);

   // check that local discrete maximum principle is satisfied at all time steps
   bool DMP_satisfied_this_step = check_max_principle(dt,false);
   DMP_satisfied_at_all_steps = DMP_satisfied_at_all_steps and DMP_satisfied_this_step;
}

/** \brief Computes min and max quantities for max principle
 */
template <int dim>
void FCT<dim>::compute_bounds(const Vector<double> &old_solution,
                              const double &dt)
{
   // initialize min and max values
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      solution_max(i) = old_solution(i);
      solution_min(i) = old_solution(i);
   }

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // find min and max values on cell
      double max_cell = old_solution(local_dof_indices[0]); // initialized to arbitrary value on cell
      double min_cell = old_solution(local_dof_indices[0]); // initialized to arbitrary value on cell
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
         double value_j = old_solution(local_dof_indices[j]);
         max_cell = std::max(max_cell, value_j);
         min_cell = std::min(min_cell, value_j);
      }

      // update the max and min values of neighborhood of each dof
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
         unsigned int i = local_dof_indices[j]; // global index
         solution_max(i) = std::max(solution_max(i), max_cell);
         solution_min(i) = std::min(solution_min(i), min_cell);
      }
   }

   // At this point, the min/max values of the old solution in the support
   // of test function i are stored in solution_min(i) and solution_max(i).
   // Now these values are multiplied by (1-dt/m(i))*sum_j(A(i,j)) and
   // added to dt/m(i)*b(i).

   // compute the upper and lower bounds for the maximum principle
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // compute sum of A_ij over row i
      // get nonzero entries of row i of A
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(low_order_ss_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);
      // add nonzero entries to get the row sum
      double row_sum = 0.0;
      for (unsigned int k = 0; k < n_col; ++k)
         row_sum += row_values[k];
      
      // compute the max and min values for the maximum principle
      solution_max(i) = solution_max(i)*(1.0 - dt/lumped_mass_matrix(i,i)*row_sum)
         + dt/lumped_mass_matrix(i,i)*ss_rhs(i);
      solution_min(i) = solution_min(i)*(1.0 - dt/lumped_mass_matrix(i,i)*row_sum)
         + dt/lumped_mass_matrix(i,i)*ss_rhs(i);
   }
}

/** \brief Computes min and max quantities for steady-state max principle
 */
template <int dim>
void FCT<dim>::compute_steady_state_max_principle_bounds()
{
   // initialize min and max values
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      solution_max(i) = new_solution(i);
      solution_min(i) = new_solution(i);
   }

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // find min and max values on cell
      double max_cell = new_solution(local_dof_indices[0]); // initialized to arbitrary value on cell
      double min_cell = new_solution(local_dof_indices[0]); // initialized to arbitrary value on cell
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
         double value_j = new_solution(local_dof_indices[j]);
         max_cell = std::max(max_cell, value_j);
         min_cell = std::min(min_cell, value_j);
      }

      // update the max and min values of neighborhood of each dof
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
         unsigned int i = local_dof_indices[j]; // global index
         solution_max(i) = std::max(solution_max(i), max_cell);
         solution_min(i) = std::min(solution_min(i), min_cell);
      }
   }

   // compute the upper and lower bounds for the maximum principle
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // compute sum of A_ij over row i
      // get nonzero entries of row i of A
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(low_order_ss_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);
      // add nonzero entries to get the row sum
      double row_sum = 0.0;
      for (unsigned int k = 0; k < n_col; ++k)
         row_sum += row_values[k];
      
      // compute the max and min values for the maximum principle
      solution_max(i) = solution_max(i)*(1.0 - row_sum/low_order_ss_matrix(i,i))
         + ss_rhs(i)/low_order_ss_matrix(i,i);
      solution_min(i) = solution_min(i)*(1.0 - row_sum/low_order_ss_matrix(i,i))
         + ss_rhs(i)/low_order_ss_matrix(i,i);
   }
}

/** \brief Assembles the high-order coefficient matrix \f$\mathcal{A}\f$.
 *  \param [in] dt current time step size
 */
template <int dim>
void FCT<dim>::compute_flux_corrections(const double &dt)
{
   // reset A to zero
   flux_correction_matrix = 0;

   // add A1 matrix to A by looping over cells
   //--------------------------------------------------------
   std::vector<unsigned int> local_dof_indices(dofs_per_cell);
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (unsigned int i_cell = 0; cell != endc; ++cell, ++i_cell)
   {
      // get global dof indices for local dofs
      cell->get_dof_indices(local_dof_indices);

      // get low-order and high-order viscosities for cell
      double nu_L = low_order_viscosity (i_cell);
      double nu_H = high_order_viscosity(i_cell);
      // compute cell volume, needed for viscous bilinear form bK(phi_j,phi_i)
      double cell_volume = cell->measure();

      for (unsigned int ii = 0; ii < dofs_per_cell; ++ii)
      {
         // get global index for local dof i
         unsigned int i = local_dof_indices[ii];
         
         for (unsigned int jj = 0; jj < dofs_per_cell; ++jj)
         {
            // get global index for local dof j
            unsigned int j = local_dof_indices[jj];

            // compute viscous bilinear form bK(phi_j,phi_i)
            double bK;
            if (ii == jj)
               bK = cell_volume;
            else
               bK = cell_volume / (1.0 - dofs_per_cell);

            // compute cell contribution to A1(i,j)
            double A1_ij = dt*(nu_L - nu_H)*(old_solution(j) - old_solution(i))*bK;
            // add it to A(i,j)
            flux_correction_matrix.add(i,j,A1_ij);
         }
      }
   }
   
   // add A2 matrix to A by looping over all degrees of freedom
   //----------------------------------------------------------
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // Note that the high-order right hand side (G) is stored in system_rhs,
      // and auxiliary_mass_matrix is the "B" matrix

      // get nonzero entries in row i of B
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(auxiliary_mass_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);

      // loop over nonzero entries in row i of B
      for (unsigned int k = 0; k < n_col; ++k)
      {
         // get the column index of nonzero entry k
         unsigned int j = row_indices[k];

         // add A2(i,j) to A(i,j)
         flux_correction_matrix.add(i,j,dt*(auxiliary_mass_matrix(j,i)*system_rhs(i)
                                                  -auxiliary_mass_matrix(i,j)*system_rhs(j)));
      }
   }
}

/** \brief Computes the limiting coefficient vectors \f$R^+\f$ and \f$R^-\f$,
 *         used for computing the high-order solution from the low-order
 *         solution.
 */
template <int dim>
void FCT<dim>::compute_limiting_coefficients()
{
   // reset flux correction vector
   flux_correction_vector = 0;

   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // for Dirichlet nodes, set R_minus and R_plus to 1 because no limiting is needed
      if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) != dirichlet_nodes.end())
      {
         R_minus(i) = 1.0;
         R_plus(i) = 1.0;
      }
      else
      {
         // get nonzero entries in row i of high-order coefficient matrix A
         std::vector<double>       row_values;
         std::vector<unsigned int> row_indices;
         unsigned int              n_col;
         get_matrix_row(flux_correction_matrix,
                        i,
                        row_values,
                        row_indices,
                        n_col);
   
         double P_plus_i  = 0.0;
         double P_minus_i = 0.0;
         // compute P_plus, P_minus
         for (unsigned int k = 0; k < n_col; ++k)
         {
            // get value of nonzero entry k
            double Fij = row_values[k];
   
            P_plus_i  += std::max(0.0,Fij);
            P_minus_i += std::min(0.0,Fij);
         }
   
         // compute R_plus(i)
         if (P_plus_i != 0.0)
            R_plus(i) = std::min(1.0, Q_plus(i)/P_plus_i);
         else
            R_plus(i) = 1.0;
   
         // compute R_minus(i)
         if (P_minus_i != 0.0)
            R_minus(i) = std::min(1.0, Q_minus(i)/P_minus_i);
         else
            R_minus(i) = 1.0;
      }
   }

   // if user chose not to limit, set R+ and R- = 1 so that limiting
   // coefficients will equal 1 and thus no limiting will occur
   if (parameters.do_not_limit)
   {
      for (unsigned int i = 0; i < n_dofs; ++i)
      {
         R_minus(i) = 1.0;
         R_plus(i)  = 1.0;
      }
   }

   // compute limited flux correction sum
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // get values and indices of nonzero entries in row i of coefficient matrix A
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(flux_correction_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);

      // perform flux correction for dof i
      // Note that flux correction sum is sum_{j in I(S_i)} of Lij*Fij.
      for (unsigned int k = 0; k < n_col; ++k) {
         unsigned int j = row_indices[k];
         double Fij     = row_values[k];
         // compute limiting coefficient Lij
         double Lij;
         if (Fij >= 0.0)
            Lij = std::min(R_plus(i),R_minus(j));
         else
            Lij = std::min(R_minus(i),R_plus(j));

         // add Lij*Fij to flux correction sum
         flux_correction_vector(i) += Lij*Fij;
      }
   }
}

/** \brief Gets the values and indices of nonzero elements in a row of a sparse matrix.
 *  \param [in] matrix sparse matrix whose row will be retrieved
 *  \param [in] i index of row to be retrieved
 *  \param [out] row_values vector of values of nonzero entries of row i
 *  \param [out] row_indices vector of indices of nonzero entries of row i
 *  \param [out] n_col number of nonzero entries of row i
 */
template <int dim>
void FCT<dim>::get_matrix_row(const SparseMatrix<double>      &matrix,
                                           const unsigned int              &i,
                                                 std::vector<double>       &row_values,
                                                 std::vector<unsigned int> &row_indices,
                                                 unsigned int              &n_col
                                          )
{
    // get first and one-past-last iterator for row
    SparseMatrix<double>::const_iterator matrix_iterator     = matrix.begin(i);
    SparseMatrix<double>::const_iterator matrix_iterator_end = matrix.end(i);

    // compute number of entries in row and then allocate memory
    n_col = matrix_iterator_end - matrix_iterator;
    row_values .resize(n_col);
    row_indices.resize(n_col);

    // loop over columns in row
    for(unsigned int k = 0; matrix_iterator != matrix_iterator_end; ++matrix_iterator, ++k)
    {
      row_values[k]  = matrix_iterator->value();
      row_indices[k] = matrix_iterator->column();
    }
}

/** \brief check that the local discrete max principle is satisfied.
 */
template<int dim>
bool FCT<dim>::check_max_principle(const double &dt,
                                                const bool   &using_high_order)
{
   // machine precision for floating point comparisons
   const double machine_tolerance = 1.0e-12;

   // now set new precision
   std::cout.precision(15);

   // check that each dof value is bounded by its neighbors
   bool local_max_principle_satisfied = true;
   for (unsigned int i = 0; i < n_dofs; ++i) {
      // check if dof is a Dirichlet node or not - if it is, don't check max principle
      if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) == dirichlet_nodes.end())
      {
         double value_i = new_solution(i);
         // check lower bound
         if (value_i < solution_min(i) - machine_tolerance) {
            local_max_principle_satisfied = false;
            // determine which condition was violated
            if (using_high_order)
            {
               std::cout << "      Max principle lower bound violated with dof "
                  << i << " of high-order solution: " << std::scientific
                  << value_i << " < " << solution_min(i) << std::endl;
            } else {
               std::cout << "      Max principle lower bound violated with dof "
                  << i << " of low-order solution: " << std::scientific
                  << value_i << " < " << solution_min(i) << std::endl;
               debug_max_principle_low_order (i,dt);
            }
         }
         // check upper bound
         if (value_i > solution_max(i) + machine_tolerance) {
            local_max_principle_satisfied = false;
            // determine which condition was violated
            if (using_high_order)
            {
               std::cout << "      Max principle upper bound violated with dof "
                  << i << " of high-order solution: " << std::scientific
                  << value_i << " > " << solution_max(i) << std::endl;
            } else {
               std::cout << "      Max principle upper bound violated with dof "
                  << i << " of low-order solution: " << std::scientific
                  << value_i << " > " << solution_max(i) << std::endl;
               debug_max_principle_low_order (i,dt);
            }
         }
      }
   }

   // restore default precision and format
   std::cout.unsetf(std::ios_base::floatfield);

   // exit if max principle was violated
   if (not local_max_principle_satisfied)
   {
      std::cout << "Program terminated due to max-principle violation." << std::endl;
      std::exit(0);
   }

   return local_max_principle_satisfied;
}

/** \brief Debugging function used for determining why the maximum
 *         principle is failing for the low-order method;
 *         examines the conditions that are required to be met for
 *         the principle to be satisfied.
 */
template <int dim>
void FCT<dim>::debug_max_principle_low_order(const unsigned int &i,
                                                          const double &dt)
{
   // flag to determine if no conditions were violated
   bool condition_violated = false;

   // get nonzero entries of row i of A
   std::vector<double>       row_values;
   std::vector<unsigned int> row_indices;
   unsigned int              n_col;
   get_matrix_row(low_order_ss_matrix,
                  i,
                  row_values,
                  row_indices,
                  n_col);

   // diagonal entry
   double Aii = low_order_ss_matrix(i,i);

   // check that system matrix diagonal elements are non-negative
   if (Aii < 0.0)
   {
      std::cout << "         DEBUG: diagonal element is negative: " << Aii << std::endl;
      condition_violated = true;
   }
   
   // loop over nonzero entries in row i to compute off-diagonal sum
   // for diagonal dominance check and also check that each off-diagonal
   // element is non-positive
   double off_diagonal_sum = 0.0;
   for (unsigned int k = 0; k < n_col; ++k)
   {
      unsigned int j = row_indices[k];
      double Aij = row_values[k];

      if (j != i)
      {
         // add to off-diagonal sum
         off_diagonal_sum += std::abs(Aij);
         // check that system matrix off-diagonal elements are non-positive
         if (Aij > 0.0)
         {
            std::cout << "         DEBUG: off-diagonal element (" << i << "," << j << ") is positive: " << Aij << std::endl;
            condition_violated = true;
         }
      }
   }

   // check that system matrix is diagonally dominant
   if (Aii < off_diagonal_sum)
   {
      std::cout << "         DEBUG: row is not diagonally dominant: Aii = "
         << Aii << ", off-diagonal sum = " << off_diagonal_sum << std::endl;
      condition_violated = true;
   }
     
   // check that CFL condition is satisfied if transient problem
   if (!parameters.is_steady_state)
   {
      double cfl = dt/lumped_mass_matrix(i,i)*low_order_ss_matrix(i,i);
      if (cfl > 1.0)
      {
         std::cout << "         DEBUG: row does not satisfy CFL condition: CFL = " << cfl << std::endl;
         condition_violated = true;
      }
   }

   // report if no conditions were violated
   if (not condition_violated)
      std::cout << "         DEBUG: No checks returned flags; deeper debugging is necessary." << std::endl;
}
