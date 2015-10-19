/**
 * \brief Constructor.
 */
template <int dim>
BoundaryConditions<dim>::BoundaryConditions(FESystem)
{
  FEValues<dim> fe_values(this->fe,
                          this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
}

/**
 * \brief Applies boundary conditions for a cell.
 *
 * \param[in] cell cell iterator
 * \param[in] fe_values_cell FE values for cell, in case cell values are needed
 * \param[in] solution solution vector
 * \param[inout] cell_residual residual vector for cell
 */
template <int dim>
void BoundaryConditions<dim>::apply(const Cell & cell,
 const FEValues<dim> & fe_values_cell,
    const Vector<double> & solution,
    Vector<double> & cell_residual)
{
  // loop over faces
  for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
  {
    // determine if face is interior
    typename DoFHandler<dim>::face_iterator face = cell->face(iface);
    if (face->at_boundary() == true)
    {
      // reinitialize FE values
      fe_values_face.reinit(cell, iface);

      // apply boundary conditions for this boundary face
      apply_boundary_condition(cell, fe_values_cell, fe_values_face, solution, cell_residual);
    }
  }
}