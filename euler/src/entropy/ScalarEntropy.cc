/**
 * \file ScalarEntropy.cc
 * \brief Provides the function definitions for the ScalarEntropy class.
 */

/**
 * \brief Constructor.
 */
template <int dim>
ScalarEntropy<dim>::ScalarEntropy()
{
}

/**
 * \brief Computes the max entropy jump in a cell.
 *
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 */
/*
template <int dim>
double ConservationLaw<dim>::compute_max_entropy_jump(
  const Vector<double> & solution, const Cell & cell) const
{
  FEFaceValues<dim> fe_values_face(*fe,
                                   *face_quadrature,
                                   update_values | update_gradients |
                                     update_normal_vectors);

  std::vector<Tensor<1, dim>> gradients_face(n_q_points_face);
  std::vector<Tensor<1, dim>> gradients_face_neighbor(n_q_points_face);
  std::vector<Tensor<1, dim>> normal_vectors(n_q_points_face);
  Vector<double> entropy(n_q_points_face);

  double max_jump_in_cell = 0.0;

  // loop over faces
  for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
  {
    // determine if face is interior
    typename DoFHandler<dim>::face_iterator face = cell->face(iface);
    if (face->at_boundary() == false)
    {
      // get gradients from this cell
      fe_values_face.reinit(cell, iface);
      fe_values_face.get_function_gradients(solution, gradients_face);

      // compute entropy at each quadrature point on face
      compute_entropy(solution, fe_values_face, entropy);

      // get normal vectors
      normal_vectors = fe_values_face.get_all_normal_vectors();

      // get iterator to neighboring cell and face index of this face
      Assert(cell->neighbor(iface).state() == IteratorState::valid,
             ExcInternalError());
      typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(iface);
      const unsigned int ineighbor = cell->neighbor_of_neighbor(iface);
      Assert(ineighbor < faces_per_cell, ExcInternalError());

      // get gradients from neighboring cell
      fe_values_face.reinit(neighbor, ineighbor);
      fe_values_face.get_function_gradients(solution, gradients_face_neighbor);

      // loop over face quadrature points to determine max jump on face
      double max_jump_on_face = 0.0;
      for (unsigned int q = 0; q < n_q_points_face; ++q)
      {
        // compute difference in gradients across face
        gradients_face[q] -= gradients_face_neighbor[q];
        double jump_on_face =
          std::abs((gradients_face[q] * normal_vectors[q]) * entropy[q]);
        max_jump_on_face = std::max(max_jump_on_face, jump_on_face);
      }

      // update max jump in cell
      max_jump_in_cell = std::max(max_jump_in_cell, max_jump_on_face);
    }
  }

  return max_jump_in_cell;
}
*/
