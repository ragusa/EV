/** \brief Constructor.
 */
template<int dim>
LowOrderViscosity<dim>::LowOrderViscosity(
   const unsigned int          n_cells,
   const unsigned int          dofs_per_cell,
   const DoFHandler<dim>      &dof_handler,
   const ConstraintMatrix     &constraints,
   const SparseMatrix<double> &inviscid_matrix,
         SparseMatrix<double> &diffusion_matrix,
         SparseMatrix<double> &total_matrix) :

   Viscosity<dim>(n_cells,dofs_per_cell,dof_handler,constraints)//,
   //dof_handler(&dof_handler)
{
   // create sparse matrix for viscous bilinear forms
   CompressedSparsityPattern c_sparsity_pattern(this->n_dofs);
   DoFTools::make_sparsity_pattern(dof_handler, c_sparsity_pattern);
   sparsity_pattern.copy_from(c_sparsity_pattern);
   viscous_bilinear_forms.reinit(sparsity_pattern);

   // compute viscous bilinear forms
   compute_viscous_bilinear_forms();

   // compute low-order viscosity
   compute_low_order_viscosity(inviscid_matrix);

   // compute diffusion matrix
   this->compute_diffusion_matrix(diffusion_matrix);

   // add diffusion matrix
   this->add_diffusion_matrix(inviscid_matrix,diffusion_matrix,total_matrix);

/*
   // allocate memory for vectors
   solution_min.reinit(dof_handler.n_dofs());
   solution_max.reinit(dof_handler.n_dofs());
*/
}

/** \brief Destructor.
 */
template<int dim>
LowOrderViscosity<dim>::~LowOrderViscosity() {
}

/** \brief Computes the low-order viscosity for each cell.
 */
template <int dim>
void LowOrderViscosity<dim>::compute_low_order_viscosity(
   const SparseMatrix<double> &inviscid_matrix)
{
   // local dof indices
   std::vector<unsigned int> local_dof_indices (this->dofs_per_cell);

   // loop over cells to compute low-order viscosity at each quadrature point
   typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler->begin_active(),
                                                  endc = this->dof_handler->end();
   unsigned int i_cell = 0;
   for (; cell != endc; ++cell, ++i_cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // initialize to zero for max() function
      double low_order_viscosity_cell = 0.0;
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
         for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
            if (i != j)
               low_order_viscosity_cell = std::max(low_order_viscosity_cell,
                  std::max(0.0,inviscid_matrix(local_dof_indices[i],local_dof_indices[j]))/
                  (-viscous_bilinear_forms(local_dof_indices[i],local_dof_indices[j])));

      // store computed viscosity value
      this->viscosity(i_cell) = low_order_viscosity_cell;
   }
}

/** \brief Computes viscous bilinear forms, to be used in the computation of
 *         maximum-principle preserving low order viscosity.
 *
 *         Each element of the resulting matrix, \f$B_{i,j}\f$ is computed as
 *         follows:
 *         \f[
 *            B_{i,j} = \sum_{K\subset S_{i,j}}b_K(\varphi_i,\varphi_j)
 *         \f]
 */
template <int dim>
void LowOrderViscosity<dim>::compute_viscous_bilinear_forms()
{
   viscous_bilinear_forms = 0; // zero out matrix

   // local dof indices
   std::vector<unsigned int> local_dof_indices (this->dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler->begin_active(),
                                                  endc = this->dof_handler->end();
   for (; cell != endc; ++cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // query cell volume
      double cell_volume = cell->measure();

      // add local bilinear forms to global matrix
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i) {
         for (unsigned int j = 0; j < this->dofs_per_cell; ++j) {
            double b_cell = 0.0;
            if (j == i) {
               b_cell = cell_volume;
            } else {
               b_cell = -1.0/(this->dofs_per_cell-1.0)*cell_volume;
            }
            viscous_bilinear_forms.add(local_dof_indices[i],
                                       local_dof_indices[j],
                                       b_cell);
         }
      }
   }
}

/** \brief Computes min and max quantities for max principle
 */
/*
template <int dim>
void LowOrderViscosity<dim>::compute_bounds(const Vector<double>       &old_solution,
                                            const SparseMatrix<double> &low_order_ss_matrix,
                                            const Vector<double>       &ss_rhs,
                                            const double               &dt)
{
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      solution_max(i) = old_solution(i);
      solution_min(i) = old_solution(i);
   }

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler->begin_active(),
                                                  endc = dof_handler->end();
   std::vector<unsigned int> local_dof_indices(dofs_per_cell);
   for (; cell != endc; ++cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // find min and max values on cell
      double max_cell = old_solution(local_dof_indices[0]);
      double min_cell = old_solution(local_dof_indices[0]);
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
      solution_max(i) = solution_max(i)*(1.0 - dt/(*lumped_mass_matrix)(i,i)*row_sum)
         + dt/(*lumped_mass_matrix)(i,i)*ss_rhs(i);
      solution_min(i) = solution_min(i)*(1.0 - dt/(*lumped_mass_matrix)(i,i)*row_sum)
         + dt/(*lumped_mass_matrix)(i,i)*ss_rhs(i);
   }
}
*/

/** \brief Outputs bounds to files.
 */
/*
template<int dim>
void LowOrderViscosity<dim>::output_bounds(const PostProcessor<dim> &postprocessor) const
{
   postprocessor.output_solution(solution_min,
                                 *dof_handler,
                                 "DMPmin");
   postprocessor.output_solution(solution_max,
                                 *dof_handler,
                                 "DMPmax");
}
*/
