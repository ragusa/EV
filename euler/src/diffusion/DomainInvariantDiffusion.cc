/**
 * \file DomainInvariantDiffusion.cc
 * \brief Provides the function definitions for the DomainInvariantDiffusion
 * class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
DomainInvariantDiffusion<dim>::DomainInvariantDiffusion()
{
}

/**
 * \brief Updates the diffusion matrix.
 *
 * \param[in] solution solution vector
 */
template <int dim>
void DomainInvariantDiffusion<dim>::update(const Vector<double> & solution)
{
  // reset diffusion matrix to zero
  this->diffusion_matrix = 0;

  // cell iterator
  Cell cell = dof_handler->begin_active(), endc = dof_handler->end();

  // reset line flags
  triangulation->clear_user_flags_line();

  // dof indices for a line
  std::vector<unsigned int> local_dof_indices(2 * n_components);

  // loop over cells
  unsigned int i_cell = 0;
  unsigned int i_line = 0;
  for (cell = dof_handler->begin_active(); cell != endc; ++cell, ++i_cell)
  {
    // loop over lines of cell
    for (unsigned int line = 0; line < GeometryInfo<dim>::lines_per_cell;
         ++line, ++i_line)
    {
      // skip line if it has already been traversed
      if (!cell->line(line)->user_flag_set())
      {
        // mark line to signal it has been traversed
        cell->line(line)->set_user_flag();

        // get dof indices on line
        cell->line(line)->get_dof_indices(local_dof_indices);
        unsigned int i = local_dof_indices[0];
        unsigned int j = local_dof_indices[1];

        // compute diffusion entry value
        double Dij = -std::max(max_speed_ij * cij, max_speed_ji * cji);

        // set diffusion matrix entries
        this->diffusion_matrix.set(i, j, Dij) this->diffusion_matrix.set(
          j, i, Dij)
      }
    }
  }
}
