/**
 * Constructor.
 *
 * @param[in] n_cells number of cells
 * @param[in] dofs_per_cell number of dofs in each cell
 * @param[in] dof_handler dof handler object
 * @param[in] constraints constraints object
 */
template <int dim>
Viscosity<dim>::Viscosity(const unsigned int n_cells,
                          const unsigned int dofs_per_cell,
                          const DoFHandler<dim> & dof_handler,
                          const ConstraintMatrix & constraints)
  : viscosity(n_cells),
    n_cells(n_cells),
    n_dofs(dof_handler.n_dofs()),
    dofs_per_cell(dofs_per_cell),
    dof_handler(&dof_handler),
    constraints(&constraints)
{
}

/**
 * Destructor.
 */
template <int dim>
Viscosity<dim>::~Viscosity()
{
}

/**
 * Adds a diffusion matrix to an inviscid matrix.
 *
 * @param[in] inviscid_matrix inviscid matrix
 * @param[in] diffusion_matrix viscous matrix
 * @param[out] total_matrix inviscid matrix plus viscous matrix
 */
template <int dim>
void Viscosity<dim>::add_diffusion_matrix(
  const SparseMatrix<double> & inviscid_matrix,
  SparseMatrix<double> & diffusion_matrix,
  SparseMatrix<double> & total_matrix)
{
  // copy inviscid component of steady-state matrix
  total_matrix.copy_from(inviscid_matrix);

  // add diffusion matrix
  total_matrix.add(1.0, diffusion_matrix);
}

/**
 * Computes a graph-theoretic diffusion matrix.
 *
 * @param[out] diffusion_matrix viscous matrix to be computed
 */
template <int dim>
void Viscosity<dim>::compute_diffusion_matrix(
  SparseMatrix<double> & diffusion_matrix)
{
  // reset diffusion matrix
  diffusion_matrix = 0;

  // cell matrix
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  // loop over cells
  unsigned int i_cell = 0;
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler->begin_active(),
                                                 endc = dof_handler->end();
  for (; cell != endc; ++cell, ++i_cell)
  {
    // reset cell matrix to zero
    cell_matrix = 0;

    // compute cell volume
    double cell_volume = cell->measure();

    // compute cell contribution to global system matrix
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        double viscous_bilinear_form;
        if (j == i)
          viscous_bilinear_form = cell_volume;
        else
          viscous_bilinear_form = cell_volume / (1.0 - dofs_per_cell);

        cell_matrix(i, j) += viscosity(i_cell) * viscous_bilinear_form;
      }

    // aggregate local matrix and rhs to global matrix and rhs
    std::vector<unsigned int> local_dof_indices(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
    constraints->distribute_local_to_global(
      cell_matrix, local_dof_indices, diffusion_matrix);
  }
}

/**
 * Gets a single viscosity value associated with a cell.
 *
 * @param[in] i cell index
 */
template <int dim>
double Viscosity<dim>::get_viscosity_value(const unsigned int i) const
{
  return viscosity(i);
}

/**
 * Outputs viscosities to file.
 *
 * @param[in] output_file name of output file
 */
template <int dim>
void Viscosity<dim>::output_viscosity(const std::string output_file) const
{
  // add viscosities to data out object
  DataOut<dim> visc_out;
  visc_out.attach_dof_handler(*dof_handler);
  visc_out.add_data_vector(viscosity, "Viscosity", DataOut<dim>::type_cell_data);

  // create output filestream
  std::ofstream viscosity_outstream(output_file.c_str());

  // build patches and write to file
  visc_out.build_patches();
  if (dim == 1)
    visc_out.write_gnuplot(viscosity_outstream);
  else
    visc_out.write_vtk(viscosity_outstream);
}
