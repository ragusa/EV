/**
 * \file DMPLowOrderViscosity.cc
 * \brief Provides the function definitions for the DMPLowOrderViscosity class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
DMPLowOrderViscosity<dim>::DMPLowOrderViscosity(
  const DoFHandler<dim> & dof_handler_,
  const SparseMatrix<double> & inviscid_matrix_,
  const unsigned int & dofs_per_cell_)
  : Viscosity<dim>(dof_handler_),
    inviscid_matrix(&inviscid_matrix_),
    dofs_per_cell(dofs_per_cell_)
{
  // reinitialize
  reinitialize();
}

/**
 * \brief Reinitializes.
 */
template <int dim>
void DMPLowOrderViscosity<dim>::reinitialize()
{
  // update number of degrees of freedom
  this->n_dofs = this->dof_handler->n_dofs();

  // create sparse matrix for viscous bilinear forms
  DynamicSparsityPattern dsp(this->n_dofs);
  DoFTools::make_sparsity_pattern(*(this->dof_handler), dsp);
  sparsity_pattern.copy_from(dsp);
  viscous_bilinear_forms.reinit(sparsity_pattern);

  // compute viscous bilinear forms
  compute_viscous_bilinear_forms();
}

/**
 * \brief Computes viscous bilinear forms, to be used in the computation of
 *        maximum-principle preserving low order viscosity.
 *
 * Each element of the resulting matrix, \f$B_{i,j}\f$ is computed as
 * follows:
 * \f[
 *   B_{i,j} = \sum_{K\subset S_{i,j}}b_K(\varphi_i,\varphi_j)
 * \f]
 */
template <int dim>
void DMPLowOrderViscosity<dim>::compute_viscous_bilinear_forms()
{
  viscous_bilinear_forms = 0; // zero out matrix

  // local dof indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells
  Cell cell = this->dof_handler->begin_active(), endc = this->dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    // query cell volume
    double cell_volume = cell->measure();

    // add local bilinear forms to global matrix
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        double b_cell = 0.0;
        if (j == i)
          b_cell = cell_volume;
        else
          b_cell = -1.0 / (dofs_per_cell - 1.0) * cell_volume;

        viscous_bilinear_forms.add(
          local_dof_indices[i], local_dof_indices[j], b_cell);
      }
    }
  }
}

/**
 * \brief Computes the maximum-principle preserving first order viscosity
 *        for each cell.
 */
template <int dim>
void DMPLowOrderViscosity<dim>::update(const Vector<double> &,
                                       const Vector<double> &,
                                       const double &,
                                       const unsigned int &)
{
  // local dof indices
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // loop over cells to compute first order viscosity at each quadrature point
  Cell cell = this->dof_handler->begin_active(), endc = this->dof_handler->end();
  for (; cell != endc; ++cell)
  {
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    this->values[cell] = 0.0;
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        if (i != j)
        {
          this->values[cell] =
            std::max(this->values[cell],
                     std::max(0.0,
                              (*inviscid_matrix)(local_dof_indices[i],
                                                 local_dof_indices[j])) /
                       (-viscous_bilinear_forms(local_dof_indices[i],
                                                local_dof_indices[j])));
        }
      }
    }
  }
}
