/**
 * \brief Constructor.
 */
template <int dim>
TransportExecutioner<dim>::TransportExecutioner(
  const TransportRunParameters<dim> & parameters_,
  TransportProblemParameters<dim> & problem_parameters_,
  Triangulation<dim> & triangulation_,
  PostProcessor<dim> & postprocessor_)
  : cout1(std::cout, parameters_.verbosity_level >= 1),
    cout2(std::cout, parameters_.verbosity_level >= 2),
    parameters(parameters_),
    problem_parameters(&problem_parameters_),
    triangulation(&triangulation_),
    fe(FE_Q<dim>(parameters.degree), 1),
    flux(0),
    dof_handler(triangulation_),
    dofs_per_cell(fe.dofs_per_cell),
    n_cells(triangulation_.n_active_cells()),
    cell_quadrature(parameters.n_quadrature_points),
    face_quadrature(parameters.n_quadrature_points),
    n_q_points_cell(cell_quadrature.size()),
    transport_direction(problem_parameters_.transport_direction),
    transport_speed(problem_parameters_.transport_speed),
    cross_section_function(&(problem_parameters_.cross_section_function)),
    source_function(&(problem_parameters_.source_function)),
    incoming_function(problem_parameters_.dirichlet_function),
    domain_volume(problem_parameters_.domain_volume),
    linear_solver(parameters.linear_solver_type,
                  constraints,
                  dof_handler,
                  incoming_function,
                  1),
    nonlinear_solver(parameters_, linear_solver, constraints),
    postprocessor(&postprocessor_)
{
  // distribute dofs
  dof_handler.distribute_dofs(fe);

  // set number of dofs
  n_dofs = dof_handler.n_dofs();

  // get number of dofs
  const unsigned int n_dofs = dof_handler.n_dofs();

  // clear constraint matrix and make hanging node constraints
  constraints.clear();
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  // create sparsity pattern
  DynamicSparsityPattern dsp(n_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  constrained_sparsity_pattern.copy_from(dsp);

  // initialize sparse matrices
  system_matrix.reinit(constrained_sparsity_pattern);
  inviscid_ss_matrix.reinit(constrained_sparsity_pattern);
  low_order_ss_matrix.reinit(constrained_sparsity_pattern);
  high_order_ss_matrix.reinit(constrained_sparsity_pattern);
  low_order_diffusion_matrix.reinit(constrained_sparsity_pattern);
  high_order_diffusion_matrix.reinit(constrained_sparsity_pattern);

  // initialize vectors
  system_rhs.reinit(n_dofs);
  ss_rhs.reinit(n_dofs);
  new_solution.reinit(n_dofs);
  cumulative_antidiffusion.reinit(n_dofs);

  // determine Dirichlet nodes
  getDirichletNodes();

  // determine low-order viscosity and artificial diffusion
  switch (parameters.low_order_scheme)
  {
    case LowOrderScheme::dmp:
      low_order_viscosity_type = ViscosityType::DMP;
      low_order_diffusion_type = DiffusionType::graphtheoretic;
      break;
    default:
      throw ExcNotImplemented();
      break;
  }

  // determine high-order viscosity and artificial diffusion
  switch (parameters.high_order_scheme)
  {
    case HighOrderScheme::galerkin:
      entropy_viscosity_type = ViscosityType::none;
      high_order_viscosity_type = ViscosityType::none;
      high_order_diffusion_type = DiffusionType::none;
      break;
    case HighOrderScheme::entropy_visc:
      entropy_viscosity_type = ViscosityType::entropy;
      high_order_viscosity_type = ViscosityType::high;
      high_order_diffusion_type = low_order_diffusion_type;
      break;
    default:
      throw ExcNotImplemented();
      break;
  }

  // set unneeded viscosities and diffusion to none
  if (parameters.scheme == Scheme::low)
  {
    entropy_viscosity_type = ViscosityType::none;
    high_order_viscosity_type = ViscosityType::none;
    high_order_diffusion_type = DiffusionType::none;
  }
  if (parameters.scheme == Scheme::high)
  {
    if (parameters.high_order_scheme == HighOrderScheme::galerkin)
    {
      low_order_viscosity_type = ViscosityType::none;
      low_order_diffusion_type = DiffusionType::none;
    }
  }

  // create low-order viscosity
  switch (low_order_viscosity_type)
  {
    // no viscosity
    case ViscosityType::none:
    {
      low_order_viscosity =
        std::make_shared<ConstantViscosity<dim>>(0.0, dof_handler);
      break;
    }
    // DMP viscosity
    case ViscosityType::DMP:
    {
      // create low-order viscosity
      low_order_viscosity = std::make_shared<DMPLowOrderViscosity<dim>>(
        dof_handler, inviscid_ss_matrix, dofs_per_cell);
      break;
    }
    default:
    {
      throw ExcNotImplemented();
      break;
    }
  }

  // create entropy viscosity
  switch (entropy_viscosity_type)
  {
    // no viscosity
    case ViscosityType::none:
    {
      entropy_viscosity =
        std::make_shared<ConstantViscosity<dim>>(0.0, dof_handler);
      break;
    }
    // entropy viscosity
    case ViscosityType::entropy:
    {
      // create entropy
      std::shared_ptr<Entropy<dim>> entropy =
        std::make_shared<TransportEntropy<dim>>(domain_volume,
                                                dof_handler,
                                                fe,
                                                cell_quadrature,
                                                face_quadrature,
                                                *problem_parameters);

      // create entropy viscosity
      entropy_viscosity = std::make_shared<EntropyViscosity<dim>>(parameters,
                                                                  entropy,
                                                                  fe,
                                                                  dof_handler,
                                                                  cell_quadrature,
                                                                  face_quadrature,
                                                                  false);
      break;
    }
    default:
    {
      throw ExcNotImplemented();
      break;
    }
  }

  // create high-order viscosity
  switch (high_order_viscosity_type)
  {
    // no viscosity
    case ViscosityType::none:
    {
      high_order_viscosity =
        std::make_shared<ConstantViscosity<dim>>(0.0, dof_handler);
      break;
    }
    // high-order viscosity
    case ViscosityType::high:
    {
      // create high-order viscosity
      high_order_viscosity = std::make_shared<HighOrderViscosity<dim>>(
        low_order_viscosity,
        entropy_viscosity,
        parameters.use_low_order_viscosity_for_first_time_step,
        dof_handler);
      break;
    }
    default:
    {
      throw ExcNotImplemented();
      break;
    }
  }

  // create low-order and high-order artificial diffusion
  low_order_diffusion = create_artificial_diffusion(low_order_diffusion_type);
  high_order_diffusion = create_artificial_diffusion(high_order_diffusion_type);
}

/**
 * \brief Creates artificial diffusion object
 *
 * \param[in] diffusion_type type of artificial diffusion
 */
template <int dim>
std::shared_ptr<ArtificialDiffusion<dim>> TransportExecutioner<
  dim>::create_artificial_diffusion(const DiffusionType & diffusion_type)
{
  std::shared_ptr<ArtificialDiffusion<dim>> artificial_diffusion;

  // create artificial diffusion
  switch (diffusion_type)
  {
    // no diffusion
    case DiffusionType::none:
    {
      artificial_diffusion = std::make_shared<NoDiffusion<dim>>();
      break;
    }
    // graph-theoretic diffusion
    case DiffusionType::graphtheoretic:
    {
      artificial_diffusion = std::make_shared<GraphTheoreticDiffusion<dim>>(
        dof_handler, dofs_per_cell, 1);
      break;
    }
    // else
    default:
    {
      throw ExcNotImplemented();
      break;
    }
  }

  return artificial_diffusion;
}

/**
 * Assembles the inviscid steady-state matrix.
 */
template <int dim>
void TransportExecutioner<dim>::assembleInviscidSteadyStateMatrix()
{
  inviscid_ss_matrix = 0;

  // FE values, for assembly terms
  FEValues<dim> fe_values(fe,
                          cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // total cross section values at each quadrature point on cell
  std::vector<double> total_cross_section_values(n_q_points_cell);

  // cell iterator
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  // loop over cells
  unsigned int i_cell = 0;
  for (cell = dof_handler.begin_active(); cell != endc; ++cell, ++i_cell)
  {
    // initialize local matrix and rhs to zero
    cell_matrix = 0;

    // reinitialize FE values
    fe_values.reinit(cell);

    // get quadrature points on cell
    std::vector<Point<dim>> points(n_q_points_cell);
    points = fe_values.get_quadrature_points();

    // get cross section values for all quadrature points
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      total_cross_section_values[q] = cross_section_function->value(points[q]);

    // compute cell contributions to global system
    // ------------------------------------------------------------------
    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    for (unsigned int q = 0; q < n_q_points_cell; ++q)
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          // store integrals of divergence and total interaction term
          // so that they may be used in computation of max-principle
          // preserving viscosity
          cell_matrix(i, j) +=
            (
              // divergence term
              fe_values[flux].value(i, q) * transport_direction *
                fe_values[flux].gradient(j, q) +
              // total interaction term
              fe_values[flux].value(i, q) * total_cross_section_values[q] *
                fe_values[flux].value(j, q)) *
            transport_speed * fe_values.JxW(q);
        } // end j
      }   // end i
    }     // end q

    // aggregate local matrix and rhs to global matrix and rhs
    constraints.distribute_local_to_global(
      cell_matrix, local_dof_indices, inviscid_ss_matrix);
  } // end cell
}

/**
 * Assembles the steady-state rhs.
 *
 * \param[out] rhs vector which to store the steady-state rhs
 * \param[in] t time at which to evaluate rhs
 */
template <int dim>
void TransportExecutioner<dim>::assembleSteadyStateRHS(Vector<double> & rhs,
                                                       const double & t)
{
  // reset steady-state rhs
  rhs = 0;

  // set the time to be used in the source function
  source_function->set_time(t);

  // FE values, for assembly terms
  FEValues<dim> fe_values(fe,
                          cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // source values at each quadrature point on cell
  std::vector<Point<dim>> points(n_q_points_cell);
  std::vector<double> source_values(n_q_points_cell);

  // cell iterator
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  // loop over cells
  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
  {
    // initialize local rhs to zero
    cell_rhs = 0;

    // reinitialize FE values
    fe_values.reinit(cell);

    // get quadrature points on cell
    points = fe_values.get_quadrature_points();

    // get total source for all quadrature points
    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      source_values[q] = source_function->value(points[q]);

    // get local dof indices
    cell->get_dof_indices(local_dof_indices);

    for (unsigned int q = 0; q < n_q_points_cell; ++q)
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        cell_rhs(i) += fe_values[flux].value(i, q) * source_values[q] *
          transport_speed * fe_values.JxW(q);

    // aggregate local matrix and rhs to global matrix and rhs
    constraints.distribute_local_to_global(cell_rhs, local_dof_indices, rhs);
  } // end cell
}

/**
 * Returns the final solution.
 */
template <int dim>
Vector<double> TransportExecutioner<dim>::getFinalSolution() const
{
  return new_solution;
}

template <int dim>
void TransportExecutioner<dim>::print_solution() const
{
  std::cout.precision(10);
  std::cout.setf(std::ios::scientific);

  for (unsigned int i = 0; i < n_dofs; ++i)
    std::cout << new_solution[i] << std::endl;
}

/**
 * Sets the boundary indicators for each boundary face.
 *
 * The Dirichlet BC is applied only to the incoming boundary, so the transport
 * direction is compared against the normal vector of the face.
 */
template <int dim>
void TransportExecutioner<dim>::setBoundaryIndicators()
{
  // reset boundary indicators to zero
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->at_boundary())
        cell->face(face)->set_boundary_id(1);

  // FE face values
  FEFaceValues<dim> fe_face_values(fe, face_quadrature, update_normal_vectors);
  // loop over cells
  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
  {
    // loop over faces of cell
    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
    {
      // if face is at boundary
      if (cell->face(face)->at_boundary())
      {
        // reinitialize FE face values
        fe_face_values.reinit(cell, face);
        // determine if the transport flux is incoming through this face;
        //  it isn't necessary to loop over all face quadrature points because
        //  the transport direction and normal vector are the same at each
        //  quadrature point; therefore, quadrature point 0 is arbitrarily chosen
        double small = -1.0e-12;
        if (fe_face_values.normal_vector(0) * transport_direction < small)
        {
          // mark boundary as incoming flux boundary: indicator 1
          cell->face(face)->set_boundary_id(0);
        }
      }
    }
  }
}

/**
 * Gets a list of dofs subject to Dirichlet boundary conditions.
 *
 * Max principle checks are not valid for Dirichlet nodes, so these nodes
 * must be excluded from limiting and DMP checks.
 */
template <int dim>
void TransportExecutioner<dim>::getDirichletNodes()
{
  // get map of Dirichlet dof indices to Dirichlet values
  std::map<unsigned int, double> boundary_values;
  VectorTools::interpolate_boundary_values(
    dof_handler, 0, ZeroFunction<dim>(), boundary_values);

  // extract dof indices from map
  dirichlet_nodes.clear();
  for (std::map<unsigned int, double>::iterator it = boundary_values.begin();
       it != boundary_values.end();
       ++it)
    dirichlet_nodes.push_back(it->first);
}

/**
 * Applies Dirichlet boundary conditions to a linear system A*x = b.
 *
 * @param [in,out] A system matrix
 * @param [in,out] b system rhs
 * @param [in,out] x system solution
 * @param [in] t  time at which Dirichlet value function is to be evaluated
 */
template <int dim>
void TransportExecutioner<dim>::applyDirichletBC(SparseMatrix<double> & A,
                                                 Vector<double> & b,
                                                 Vector<double> & x,
                                                 const double & t)
{
  // set time for dirichlet function
  incoming_function->set_time(t);

  // create map of dofs to boundary values to be imposed
  std::map<unsigned int, double> boundary_values;
  VectorTools::interpolate_boundary_values(
    dof_handler, 0, *incoming_function, boundary_values);

  // apply boundary values to system
  MatrixTools::apply_boundary_values(boundary_values, A, x, b);
}
