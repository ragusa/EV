/**
 * \file BoundaryConditions.cc
 * \brief Provides the function definitions for the BoundaryConditions class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] fe_ finite element system
 * \param[in] face_quadrature_ quadrature for face
 */
template <int dim>
BoundaryConditions<dim>::BoundaryConditions(
  const FESystem<dim> & fe_, const QGauss<dim - 1> & face_quadrature_)
  : fe(fe_),
    face_quadrature(face_quadrature_),
    fe_values_face(fe,
                   face_quadrature,
                   update_values | update_gradients | update_normal_vectors |
                     update_quadrature_points | update_JxW_values),
    faces_per_cell(GeometryInfo<dim>::faces_per_cell),
    dofs_per_cell(fe.dofs_per_cell),
    n_quadrature_points_face(face_quadrature.size())
{
}

/**
 * \brief Loops over faces in a cell to apply its boundary conditions, if any.
 *
 * This version is for building a system \e vector.
 *
 * \pre It is assumed that the cell FE values have already been reinitialized to
 *      the cell.
 *
 * \param[in] cell cell iterator
 * \param[in] fe_values_cell FE values for cell, in case cell values are needed
 * \param[in] solution solution vector
 * \param[in] dt time step size \f$\Delta t\f$
 * \param[inout] cell_residual residual vector for cell
 */
template <int dim>
void BoundaryConditions<dim>::apply(const Cell & cell,
                                    const FEValues<dim> & fe_values_cell,
                                    const Vector<double> & solution,
                                    const double & dt,
                                    Vector<double> & cell_residual)
{
  // loop over faces
  for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
  {
    // determine if face is interior
    Face face = cell->face(iface);
    if (face->at_boundary() == true)
    {
      // reinitialize FE values
      fe_values_face.reinit(cell, iface);

      // apply boundary conditions for this boundary face
      apply_boundary_condition(
        cell, fe_values_cell, fe_values_face, solution, dt, cell_residual);
    }
  }
}

/**
 * \brief Loops over faces in a cell to apply its boundary conditions, if any.
 *
 * This version is for building a system \e matrix.
 *
 * \pre It is assumed that the cell FE values have already been reinitialized to
 *      the cell.
 *
 * \param[in] cell cell iterator
 * \param[in] fe_values_cell FE values for cell, in case cell values are needed
 * \param[in] solution solution vector
 * \param[in] dt time step size \f$\Delta t\f$
 * \param[inout] cell_matrix  matrix for cell
 */
template <int dim>
void BoundaryConditions<dim>::apply(const Cell & cell,
                                    const FEValues<dim> & fe_values_cell,
                                    const Vector<double> & solution,
                                    const double & dt,
                                    FullMatrix<double> & cell_matrix)
{
  // loop over faces
  for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
  {
    // determine if face is interior
    Face face = cell->face(iface);
    if (face->at_boundary() == true)
    {
      // reinitialize FE values
      fe_values_face.reinit(cell, iface);

      // apply boundary conditions for this boundary face
      apply_boundary_condition_matrix(
        cell, fe_values_cell, fe_values_face, solution, dt, cell_matrix);
    }
  }
}

/**
 * \brief Applies boundary condition for a face.
 *
 * \param[in] cell cell iterator
 * \param[in] fe_values_cell FE values for cell
 * \param[in] fe_values_face FE values for face
 * \param[in] solution solution vector
 * \param[in] dt time step size \f$\Delta t\f$
 * \param[inout] cell_matrix steady-state matrix for cell
 */
template <int dim>
void BoundaryConditions<dim>::apply_boundary_condition_matrix(
  const Cell &,
  const FEValues<dim> &,
  const FEFaceValues<dim> &,
  const Vector<double> &,
  const double &,
  FullMatrix<double> &)
{
  // throw exception if derived class tries to call this method but did
  // not implement it
  AssertThrow(false, ExcNotImplemented());
}
