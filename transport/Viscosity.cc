/* \brief Constructor.
 */
template <int dim>
Viscosity<dim>::Viscosity(const unsigned int      n_cells,
                          const unsigned int      dofs_per_cell,
                          const DoFHandler<dim>  &dof_handler) :
   diffusion_matrix(),
   viscosity(n_cells),
   n_cells(n_cells),
   dofs_per_cell(dofs_per_cell),
   dof_handler(&dof_handler)
{
   // clear constraint matrix and make hanging node constraints
   constraints.clear();
   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
   constraints.close();

   // create sparsity pattern for diffusion matrix
   CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern(dof_handler,
                                   compressed_sparsity_pattern,
                                   constraints,
                                   false);
   sparsity_pattern.copy_from(compressed_sparsity_pattern);

   // initialize diffusion matrix
   diffusion_matrix.reinit(sparsity_pattern);
}

/* \brief Destructor.
 */
template <int dim>
Viscosity<dim>::~Viscosity()
{
}

/* \brief Computes a steady-state matrix with diffusion.
 */
template <int dim>
void Viscosity<dim>::compute_ss_matrix(const SparseMatrix<double> &inviscid_ss_matrix,
                                             SparseMatrix<double> &ss_matrix)
{
   // copy inviscid component of steady-state matrix
   ss_matrix.copy_from(inviscid_ss_matrix);

   // compute diffusion matrix
   compute_diffusion_matrix(diffusion_matrix);

   // add diffusion matrix
   ss_matrix.add(1.0, diffusion_matrix);
}

/** \brief Computes a graph-Laplacian diffusion matrix.
 */
template <int dim>
void Viscosity<dim>::compute_diffusion_matrix()
{
   // reset diffusion matrix
   diffusion_matrix = 0;

   // cell matrix
   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

   // loop over cells
   unsigned int i_cell = 0;
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler->begin_active(),
                                                  endc = dof_handler->end();
   for (; cell != endc; ++cell, ++i_cell) {
      // reset cell matrix to zero
      cell_matrix = 0;

      // compute cell volume
      double cell_volume = cell->measure();

      // compute cell contribution to global system matrix
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
         for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            double viscous_bilinear_form;
            if (j == i)
               viscous_bilinear_form = cell_volume;
            else
               viscous_bilinear_form = cell_volume/(1.0 - dofs_per_cell);
 
            cell_matrix(i,j) += viscosity(i_cell) * viscous_bilinear_form;
         }

      // aggregate local matrix and rhs to global matrix and rhs
      std::vector<unsigned int> local_dof_indices(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix,
                                             local_dof_indices,
                                             diffusion_matrix);

   }
}
