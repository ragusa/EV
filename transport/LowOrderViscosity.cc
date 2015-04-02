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

   Viscosity<dim>(n_cells,dofs_per_cell,dof_handler,constraints)
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
