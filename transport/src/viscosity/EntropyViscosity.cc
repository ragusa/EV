/**
 * Constructor.
 *
 * @param[in] fe finite element object
 *
 */
template <int dim>
EntropyViscosity<dim>::EntropyViscosity(
  const FESystem<dim> & fe,
  const unsigned int & n_cells,
  const DoFHandler<dim> & dof_handler,
  const ConstraintMatrix & constraints,
  const QGauss<dim> & cell_quadrature,
  const QGauss<dim - 1> & face_quadrature,
  const Tensor<1, dim> & transport_direction,
  const double & transport_speed,
  const FunctionParser<dim> & cross_section_function,
  FunctionParser<dim> & source_function,
  const std::string & entropy_string,
  const std::string & entropy_derivative_string,
  const double & entropy_residual_coefficient,
  const double & jump_coefficient,
  const double & domain_volume,
  const EntropyTemporalDiscretization & temporal_discretization,
  const LowOrderViscosity<dim> & low_order_viscosity,
  const SparseMatrix<double> & inviscid_matrix,
  SparseMatrix<double> & high_order_diffusion_matrix,
  SparseMatrix<double> & total_matrix)
  : Viscosity<dim>(n_cells, fe.dofs_per_cell, dof_handler, constraints),
    fe(&fe),
    flux(0),
    n_dofs(dof_handler.n_dofs()),
    faces_per_cell(GeometryInfo<dim>::faces_per_cell),
    cell_quadrature(cell_quadrature),
    face_quadrature(face_quadrature),
    n_q_points_cell(cell_quadrature.size()),
    n_q_points_face(face_quadrature.size()),
    transport_direction(transport_direction),
    transport_speed(transport_speed),
    cross_section_function(&cross_section_function),
    source_function(&source_function),
    entropy_string(entropy_string),
    entropy_derivative_string(entropy_derivative_string),
    entropy_residual_coefficient(entropy_residual_coefficient),
    jump_coefficient(jump_coefficient),
    domain_volume(domain_volume),
    entropy_viscosity(n_cells),
    temporal_discretization(temporal_discretization),
    low_order_viscosity(&low_order_viscosity),
    inviscid_matrix(&inviscid_matrix),
    high_order_diffusion_matrix(&high_order_diffusion_matrix),
    total_matrix(&total_matrix)
{
  // initialize entropy function
  std::map<std::string, double> constants;
  entropy_function.initialize("u", entropy_string, constants, false);
  entropy_derivative_function.initialize(
    "u", entropy_derivative_string, constants, false);
}

/**
 * \brief Recomputes the high-order diffusion matrix.
 *
 * \param[in] old_solution  old solution
 * \param[in] older_solution  older solution
 * \param[in] oldest_solution  oldest solution
 * \param[in] old_dt  old time step size
 * \param[in] older_dt  older time step size
 * \param[in] time  time at which to evaluate entropy residual
 * \param[out] diffusion_matrix  diffusion matrix
 */
template <int dim>
void EntropyViscosity<dim>::recompute_high_order_diffusion_matrix(
  const Vector<double> & old_solution,
  const Vector<double> & older_solution,
  const Vector<double> & oldest_solution,
  const double & old_dt,
  const double & older_dt,
  const double & time,
  SparseMatrix<double> & diffusion_matrix)
{
  // recompute entropy viscosity
  compute_entropy_viscosity(
    old_solution, older_solution, oldest_solution, old_dt, older_dt, time);

  // compute diffusion matrix
  this->compute_diffusion_matrix(diffusion_matrix);
}

/**
 * \brief Recomputes the high-order steady-state matrix.
 *
 * \param[in] old_solution  old solution
 * \param[in] older_solution  older solution
 * \param[in] oldest_solution  oldest solution
 * \param[in] old_dt  old time step size
 * \param[in] older_dt  older time step size
 * \param[in] time  time at which to evaluate entropy residual
 * \param[out] diffusion_matrix  diffusion matrix
 * \param[out] ss_matrix  steady-state matrix
 */
template <int dim>
void EntropyViscosity<dim>::recompute_high_order_ss_matrix(
  const Vector<double> & old_solution,
  const Vector<double> & older_solution,
  const Vector<double> & oldest_solution,
  const double & old_dt,
  const double & older_dt,
  const double & time,
  SparseMatrix<double> & diffusion_matrix,
  SparseMatrix<double> & ss_matrix)
{
  // compute diffusion matrix
  recompute_high_order_diffusion_matrix(old_solution,
                                        older_solution,
                                        oldest_solution,
                                        old_dt,
                                        older_dt,
                                        time,
                                        diffusion_matrix);

  // add diffusion matrix
  this->add_diffusion_matrix(*inviscid_matrix, diffusion_matrix, ss_matrix);
}

/**
 * \brief Recomputes the high-order steady-state matrix.
 *
 * This version of the function is for steady-state computations.
 * It employs the same function used to compute entropy viscosity
 * for transient computations but supplies dummy values to time
 * step size and time parameters to trick the function into computing
 * the steady-state entropy residual.
 *
 * \param[in] solution  solution
 */
template <int dim>
void EntropyViscosity<dim>::recompute_high_order_ss_matrix(
  const Vector<double> & solution)
{
  // recompute entropy viscosity
  compute_entropy_viscosity(solution,
                            solution,
                            solution,
                            1.0,  // arbitrary value for time step size
                            1.0,  // arbitrary value for time step size
                            0.0); // arbitrary value for time

  // compute diffusion matrix
  this->compute_diffusion_matrix(*high_order_diffusion_matrix);

  // add diffusion matrix
  this->add_diffusion_matrix(
    *inviscid_matrix, *high_order_diffusion_matrix, *total_matrix);
}

/**
 * Computes the domain-averaged entropy and the max entropy deviation in the
 * domain.
 *
 * @param[in] solution solution with which to calculate entropy for the
 *            normalization constant
 */
template <int dim>
void EntropyViscosity<dim>::compute_normalization_constant(
  const Vector<double> & solution)
{
  FEValues<dim> fe_values(
    *fe, cell_quadrature, update_values | update_JxW_values);
  std::vector<double> solution_values(n_q_points_cell);
  std::vector<double> entropy_values(n_q_points_cell);
  std::vector<Point<1>> solution_local_points(n_q_points_cell);

  // compute domain-averaged entropy
  //--------------------------------
  double domain_integral_entropy = 0.0;

  // loop over cells
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          ->begin_active(),
                                                 endc = this->dof_handler->end();
  for (cell = this->dof_handler->begin_active(); cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);
    // get solution values
    fe_values[flux].get_function_values(solution, solution_values);
    // compute entropy_values values
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      solution_local_points[q][0] = solution_values[q];
    entropy_function.value_list(solution_local_points, entropy_values);
    // loop over quadrature points
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      // add contribution of quadrature point to entropy integral
      domain_integral_entropy += entropy_values[q] * fe_values.JxW(q);
    }
  }
  // domain-averaged entropy_values
  domain_averaged_entropy = domain_integral_entropy / domain_volume;

  // compute max deviation of entropy_values from domain-averaged entropy
  //--------------------------------------------------------------
  normalization_constant = 0.0;

  // loop over cells
  for (cell = this->dof_handler->begin_active(); cell != endc; ++cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);
    // get old values
    std::vector<double> old_solution_values(n_q_points_cell);
    fe_values[flux].get_function_values(solution, old_solution_values);
    // compute entropy_values values
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      solution_local_points[q][0] = old_solution_values[q];
    entropy_function.value_list(solution_local_points, entropy_values);
    // loop over quadrature points
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      // add contribution of quadrature point to entropy_values integral
      normalization_constant =
        std::max(normalization_constant,
                 std::abs(entropy_values[q] - domain_averaged_entropy));
    }
  }

  // guard against division by zero
  if (normalization_constant == 0.0)
    normalization_constant = 1.0;
}

/**
 * Computes the entropy viscosity for each cell.
 */
template <int dim>
void EntropyViscosity<dim>::compute_entropy_viscosity(
  const Vector<double> & old_solution,
  const Vector<double> & older_solution,
  const Vector<double> & oldest_solution,
  const double & old_dt,
  const double & older_dt,
  const double & time)
{
  // compute the temporal discretization constants used in entropy residual
  compute_temporal_discretization_constants(old_dt, older_dt);

  // set the time of the source function
  source_function->set_time(time);

  // compute max entropy deviation in domain
  compute_normalization_constant(old_solution);

  // cell values
  std::vector<Point<dim>> points(n_q_points_cell);
  std::vector<double> sigma(n_q_points_cell);
  std::vector<double> source(n_q_points_cell);
  std::vector<double> u_old(n_q_points_cell);
  std::vector<double> u_older(n_q_points_cell);
  std::vector<double> u_oldest(n_q_points_cell);
  std::vector<Point<1>> u_old_points(n_q_points_cell);
  std::vector<Point<1>> u_older_points(n_q_points_cell);
  std::vector<Point<1>> u_oldest_points(n_q_points_cell);
  std::vector<Tensor<1, dim>> dudx_old(n_q_points_cell);
  std::vector<Tensor<1, dim>> dudx_older(n_q_points_cell);
  std::vector<double> s_old(n_q_points_cell);
  std::vector<double> s_older(n_q_points_cell);
  std::vector<double> s_oldest(n_q_points_cell);
  std::vector<double> dsdu_old(n_q_points_cell);
  std::vector<double> dsdu_older(n_q_points_cell);

  // face values
  std::vector<Tensor<1, dim>> normal(n_q_points_face);
  std::vector<double> u_old_face(n_q_points_face);
  std::vector<double> u_old_face_neighbor(n_q_points_face);
  std::vector<Tensor<1, dim>> dudx_old_face(n_q_points_face);
  std::vector<Tensor<1, dim>> dudx_old_face_neighbor(n_q_points_face);
  std::vector<Point<1>> u_old_face_points(n_q_points_face);
  std::vector<Point<1>> u_old_face_points_neighbor(n_q_points_face);
  std::vector<double> dsdu_old_face(n_q_points_face);
  std::vector<double> dsdu_old_face_neighbor(n_q_points_face);

  // FE cell values for computing entropy
  FEValues<dim> fe_values(*fe,
                          cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  // FE face values for computing entropy jumps
  FEFaceValues<dim> fe_values_face(*fe,
                                   face_quadrature,
                                   update_values | update_gradients |
                                     update_JxW_values | update_normal_vectors);
  FEFaceValues<dim> fe_values_face_neighbor(
    *fe,
    face_quadrature,
    update_values | update_gradients | update_JxW_values | update_normal_vectors);

  // cell iterator
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          ->begin_active(),
                                                 endc = this->dof_handler->end();
  // loop over cells
  unsigned int i_cell = 0;
  for (cell = this->dof_handler->begin_active(); cell != endc; ++cell, ++i_cell)
  {
    // reinitialize FE values
    fe_values.reinit(cell);

    // get solution values and gradients
    fe_values[flux].get_function_values(old_solution, u_old);
    fe_values[flux].get_function_values(older_solution, u_older);
    fe_values[flux].get_function_values(oldest_solution, u_oldest);
    fe_values[flux].get_function_gradients(old_solution, dudx_old);
    fe_values[flux].get_function_gradients(older_solution, dudx_older);

    // get cross section and source values for all quadrature points
    points = fe_values.get_quadrature_points();
    source_function->value_list(points, source);
    cross_section_function->value_list(points, sigma);

    // compute max entropy residual in cell
    //----------------------------------------------------------------------------
    // compute entropy values at each quadrature point on cell
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      u_old_points[q][0] = u_old[q];
      u_older_points[q][0] = u_older[q];
      u_oldest_points[q][0] = u_oldest[q];
    }
    entropy_function.value_list(u_old_points, s_old);
    entropy_function.value_list(u_older_points, s_older);
    entropy_function.value_list(u_oldest_points, s_oldest);
    entropy_derivative_function.value_list(u_old_points, dsdu_old);
    entropy_derivative_function.value_list(u_older_points, dsdu_older);

    // compute entropy residual values at each quadrature point on cell
    std::vector<double> entropy_residual_values(n_q_points_cell, 0.0);
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      entropy_residual_values[q] = a_old * s_old[q] + a_older * s_older[q] +
        a_oldest * s_oldest[q] +
        b_old * (dsdu_old[q] * (transport_speed * transport_direction * dudx_old[q] +
                                sigma[q] * u_old[q] - source[q])) +
        b_older * (dsdu_older[q] * (transport_speed * transport_direction * dudx_older[q] +
                                    transport_speed * sigma[q] * u_older[q] - transport_speed * source[q]));

    // determine maximum entropy residual in cell
    double max_entropy_residual = 0.0;
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      max_entropy_residual =
        std::max(max_entropy_residual, std::max(0.0, entropy_residual_values[q]));
    }

    // compute max jump in cell
    //----------------------------------------------------------------------------
    double max_jump_in_cell = 0.0;
    double max_jump_on_face = 0.0;
    for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
    {
      typename DoFHandler<dim>::face_iterator face = cell->face(iface);
      if (face->at_boundary() == false)
      {
        Assert(cell->neighbor(iface).state() == IteratorState::valid,
               ExcInternalError());
        typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(iface);
        const unsigned int ineighbor = cell->neighbor_of_neighbor(iface);
        Assert(ineighbor < faces_per_cell, ExcInternalError());

        fe_values_face.reinit(cell, iface);
        fe_values_face_neighbor.reinit(neighbor, ineighbor);

        // get solution and gradients on face
        fe_values_face.get_function_values(old_solution, u_old_face);
        fe_values_face_neighbor.get_function_values(old_solution,
                                                    u_old_face_neighbor);
        fe_values_face.get_function_gradients(old_solution, dudx_old_face);
        fe_values_face_neighbor.get_function_gradients(old_solution,
                                                       dudx_old_face_neighbor);

        // get normal vectors
        normal = fe_values_face.get_all_normal_vectors();

        // compute entropy at each quadrature point on face
        for (unsigned int q = 0; q < n_q_points_face; ++q)
        {
          u_old_face_points[q][0] = u_old_face[q];
          u_old_face_points_neighbor[q][0] = u_old_face_neighbor[q];
        }
        entropy_derivative_function.value_list(u_old_face_points, dsdu_old_face);
        entropy_derivative_function.value_list(u_old_face_points_neighbor,
                                               dsdu_old_face_neighbor);

        // compute max jump on face
        max_jump_on_face = 0.0;
        for (unsigned int q = 0; q < n_q_points_face; ++q)
        {
          double jump_dsdn =
            normal[q] * (dsdu_old_face[q] * dudx_old_face[q] -
                         dsdu_old_face_neighbor[q] * dudx_old_face_neighbor[q]);
          double jump_on_face =
            std::abs(transport_speed * transport_direction * normal[q] * jump_dsdn);
          max_jump_on_face = std::max(max_jump_on_face, jump_on_face);
        }
      }
      max_jump_in_cell = std::max(max_jump_in_cell, max_jump_on_face);
    } // end face loop

    // compute entropy viscosity in cell
    //----------------------------------------------------------------------------
    double entropy_viscosity_cell =
      (entropy_residual_coefficient * max_entropy_residual +
       jump_coefficient * max_jump_in_cell) /
      normalization_constant;

    // get low-order viscosity for cell
    double low_order_viscosity_cell =
      low_order_viscosity->get_viscosity_value(i_cell);

    // compute high-order viscosity
    this->viscosity(i_cell) =
      std::min(entropy_viscosity_cell, low_order_viscosity_cell);
  }
}

/**
 * \brief Computes the temporal discretization constants used in the entropy
 *        residual.
 */
template <int dim>
void EntropyViscosity<dim>::compute_temporal_discretization_constants(
  const double old_dt, const double older_dt)
{
  switch (temporal_discretization)
  {
    case EntropyTemporalDiscretization::BE:
    {
      a_old = 1.0 / old_dt;
      a_older = -1.0 / old_dt;
      a_oldest = 0.0;
      b_old = 1.0;
      b_older = 0.0;
      break;
    }
    case EntropyTemporalDiscretization::CN:
    {
      a_old = 1.0 / old_dt;
      a_older = -1.0 / old_dt;
      a_oldest = 0.0;
      b_old = 0.5;
      b_older = 0.5;
      break;
    }
    case EntropyTemporalDiscretization::BDF2:
    {
      a_old = (older_dt + 2 * old_dt) / (old_dt * (older_dt + old_dt));
      a_older = -(older_dt + old_dt) / (older_dt * old_dt);
      a_oldest = old_dt / (older_dt * (older_dt + old_dt));
      b_old = 1.0;
      b_older = 0.0;
      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
    }
  }
}
