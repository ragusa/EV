/**
 * \file TransportEntropy.cc
 * \brief Provides the function definitions for the TransportEntropy class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] domain_volume_ volume of domain
 * \param[in] dof_handler_ degree of freedom handler
 * \param[in] fe_ finite element system
 * \param[in] cell_quadrature_ cell quadrature
 * \param[in] face_quadrature_ face quadrature
 * \param[in] problem_parameters_ problem parameters
 */
template <int dim>
TransportEntropy<dim>::TransportEntropy(
  const double & domain_volume_,
  const DoFHandler<dim> & dof_handler_,
  const FESystem<dim> & fe_,
  const QGauss<dim> & cell_quadrature_,
  const QGauss<dim - 1> & face_quadrature_,
  const TransportProblemParameters<dim> & problem_parameters_)
  : Entropy<dim>(true,
                 domain_volume_,
                 dof_handler_,
                 fe_,
                 cell_quadrature_,
                 face_quadrature_),
    transport_speed(problem_parameters_.transport_speed),
    transport_direction(problem_parameters_.transport_direction),
    cross_section_function(&problem_parameters_.cross_section_function),
    source_function(&problem_parameters_.source_function)
{
}

/**
 * \brief Computes the entropy residual at each quadrature point in a cell.
 *
 * \param[in] new_solution new solution
 * \param[in] old_solution old solution
 * \param[in] dt time step size
 * \param[in] cell cell iterator
 *
 * \return vector of entropy residual for each quadrature point in cell
 */
template <int dim>
std::vector<double> TransportEntropy<dim>::compute_entropy_residual(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const double & dt,
  const Cell & cell)
{
  // FE values
  FEValues<dim> fe_values(*this->fe,
                          *this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points);
  fe_values.reinit(cell);

  // get solution at quadrature points
  std::vector<double> solution_local_new(this->n_q_points_cell);
  std::vector<double> solution_local_old(this->n_q_points_cell);
  fe_values.get_function_values(new_solution, solution_local_new);
  fe_values.get_function_values(old_solution, solution_local_old);

  // compute solution gradient
  std::vector<Tensor<1, dim>> solution_gradients(this->n_q_points_cell);
  fe_values.get_function_gradients(new_solution, solution_gradients);

  // compute entropy of current and old solutions
  std::vector<double> entropy_new = compute_entropy_scalar(solution_local_new);
  std::vector<double> entropy_old = compute_entropy_scalar(solution_local_old);

  // compute entropy derivative
  std::vector<double> entropy_deriv =
    compute_entropy_derivative(solution_local_new);

  // get quadrature points on cell
  std::vector<Point<dim>> points(this->n_q_points_cell);
  points = fe_values.get_quadrature_points();

  // compute entropy residual at each quadrature point on cell
  std::vector<double> entropy_residual(this->n_q_points_cell);
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
  {
    // compute cross section and source values
    const double sigma_q = cross_section_function->value(points[q]);
    const double source_q = source_function->value(points[q]);

    // compute entropy residual
    entropy_residual[q] = (entropy_new[q] - entropy_old[q]) / dt +
      entropy_deriv[q] * transport_speed *
        (transport_direction * solution_gradients[q] +
         sigma_q * solution_local_new[q] - source_q);
  }

  return entropy_residual;
}

/**
 * \brief Computes the max entropy jump in a cell.
 *
 * The max entropy jump on a cell is computed as
 * \f[
 *   J_K^n = \max\limits_{F\in\mathcal{F}(K)} J_F^n \,,
 * \f]
 * where
 * \f[
 *   J_F^n = \|v\mathbf{\Omega}\cdot\mathbf{n}_F
 *     [[\eta'(u^n)\nabla u^n\cdot\mathbf{n}_F]]\|_{L^\infty(F)} \,.
 * \f]
 *
 * \param[in] solution  solution vector
 * \param[in] cell      cell iterator
 *
 * \return max entropy jump in cell
 */
template <int dim>
double TransportEntropy<dim>::compute_max_entropy_jump(
  const Vector<double> & solution, const Cell & cell)
{
  // FE face values
  FEFaceValues<dim> fe_face_values(*this->fe,
                                   *this->face_quadrature,
                                   update_values | update_gradients |
                                     update_normal_vectors);
  FEFaceValues<dim> fe_face_values_neighbor(
    *this->fe, *this->face_quadrature, update_values | update_gradients);

  // local solution, solution gradients, and entropy derivative
  std::vector<double> solution_local(this->n_q_points_face);
  std::vector<double> solution_local_neighbor(this->n_q_points_face);
  std::vector<Tensor<1, dim>> solution_gradients(this->n_q_points_face);
  std::vector<Tensor<1, dim>> solution_gradients_neighbor(this->n_q_points_face);
  std::vector<double> entropy_deriv(this->n_q_points_face);
  std::vector<double> entropy_deriv_neighbor(this->n_q_points_face);
  std::vector<Tensor<1, dim>> normal(this->n_q_points_face);

  // initialize max entropy jump in cell to zero
  double max_jump_in_cell = 0.0;

  // loop over faces
  for (unsigned int iface = 0; iface < this->faces_per_cell; ++iface)
  {
    // determine if face is interior
    typename DoFHandler<dim>::face_iterator face = cell->face(iface);
    if (face->at_boundary() == false)
    {
      // get iterator to neighboring cell and face index of this face
      Assert(cell->neighbor(iface).state() == IteratorState::valid,
             ExcInternalError());
      Cell neighbor = cell->neighbor(iface);
      const unsigned int iface_neighbor = cell->neighbor_of_neighbor(iface);
      Assert(iface_neighbor < this->faces_per_cell, ExcInternalError());

      // reinitialize FE face values
      fe_face_values.reinit(cell, iface);
      fe_face_values_neighbor.reinit(neighbor, iface_neighbor);

      // get solution and gradients on face
      fe_face_values.get_function_values(solution, solution_local);
      fe_face_values_neighbor.get_function_values(solution,
                                                  solution_local_neighbor);
      fe_face_values.get_function_gradients(solution, solution_gradients);
      fe_face_values_neighbor.get_function_gradients(solution,
                                                     solution_gradients_neighbor);

      // compute entropy derivative for this cell and neighbor
      std::vector<double> entropy_deriv =
        compute_entropy_derivative(solution_local);
      std::vector<double> entropy_deriv_neighbor =
        compute_entropy_derivative(solution_local_neighbor);

      // get normal vectors
      normal = fe_face_values.get_all_normal_vectors();

      // loop over face quadrature points to determine max jump on face
      double max_jump_on_face = 0.0;
      for (unsigned int q = 0; q < this->n_q_points_face; ++q)
      {
        // compute jump in normal component of entropy gradient
        const double jump_normal_entropy_gradient = normal[q] *
          (entropy_deriv[q] * solution_gradients[q] -
           entropy_deriv_neighbor[q] * solution_gradients_neighbor[q]);

        // finish computing entropy jump for face
        const double jump_on_face =
          std::abs(transport_speed * transport_direction * normal[q] *
                   jump_normal_entropy_gradient);

        // update max jump on face
        max_jump_on_face = std::max(max_jump_on_face, jump_on_face);
      }

      // update max jump in cell
      max_jump_in_cell = std::max(max_jump_in_cell, max_jump_on_face);
    }
  }

  return max_jump_in_cell;
}

/**
 * \brief Computes entropy \f$\eta\f$ at each quadrature point on cell or face.
 *
 * For now, the entropy function is assumed to be \f$\eta(u) =
 * \frac{1}{2}u^2\f$.
 *
 * \param[in] solution   solution vector
 * \param[in] fe_values  FE values
 *
 * \return vector of entropy values at each quadrature point
 */
template <int dim>
std::vector<double> TransportEntropy<dim>::compute_entropy(
  const Vector<double> & solution, const FEValuesBase<dim> & fe_values) const
{
  // get local solution
  std::vector<double> solution_local(this->n_q_points_cell);
  fe_values.get_function_values(solution, solution_local);

  // call scalar entropy function
  std::vector<double> entropy = compute_entropy_scalar(solution_local);

  return entropy;
}

/**
 * \brief Computes entropy \f$\eta\f$ at each quadrature point on cell or face.
 *
 * For now, the entropy function is assumed to be \f$\eta(u) =
 * \frac{1}{2}u^2\f$.
 *
 * \param[in] solution_local  vector of solution values
 *
 * \return vector of entropy values at each quadrature point
 */
template <int dim>
std::vector<double> TransportEntropy<dim>::compute_entropy_scalar(
  const std::vector<double> & solution_local) const noexcept
{
  // get number of quadrature points
  const unsigned int n = solution_local.size();

  // compute entropy
  std::vector<double> entropy(n);
  for (unsigned int q = 0; q < n; ++q)
    entropy[q] = 0.5 * solution_local[q] * solution_local[q];

  return entropy;
}

/**
 * \brief Computes entropy derivative \f$\eta'(u)\f$ at each quadrature point on
 *        cell or face.
 *
 * For now, the entropy function is assumed to be \f$\eta(u) =
 * \frac{1}{2}u^2\f$,
 * so the entropy derivative is \f$\eta'(u) = u\f$.
 *
 * \param[in] solution_local  vector of solution values
 *
 * \return vector of entropy derivative values at each quadrature point
 */
template <int dim>
std::vector<double> TransportEntropy<dim>::compute_entropy_derivative(
  const std::vector<double> & solution_local) const noexcept
{
  // get size of vector
  const unsigned int n = solution_local.size();

  // compute entropy
  std::vector<double> entropy_derivative(n);
  for (unsigned int q = 0; q < n; ++q)
    entropy_derivative[q] = solution_local[q];

  return entropy_derivative;
}
