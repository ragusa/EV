/**
 * Constructor.
 */
template<int dim>
Executioner<dim>::Executioner(const TransportParameters<dim> & parameters_,
  Triangulation<dim> & triangulation_,
  const Tensor<1, dim> & transport_direction_,
  const FunctionParser<dim> & cross_section_function_,
  FunctionParser<dim> & source_function_, Function<dim> & incoming_function_,
  const double & domain_volume_,
  PostProcessor<dim> & postprocessor_) :
    parameters(parameters_),
    triangulation(& triangulation_),
    fe(FE_Q<dim>(parameters.degree), 1),
    flux(0),
    dof_handler(triangulation_),
    dofs_per_cell(fe.dofs_per_cell),
    n_cells(triangulation_.n_active_cells()),
    cell_quadrature(parameters.n_quadrature_points),
    face_quadrature(parameters.n_quadrature_points),
    n_q_points_cell(cell_quadrature.size()),
    linear_solver(parameters.linear_solver_option, constraints, dof_handler,
      incoming_function_),
    transport_direction(transport_direction_),
    cross_section_function(& cross_section_function_),
    source_function(& source_function_),
    incoming_function(& incoming_function_),
    domain_volume(domain_volume_),
    postprocessor(& postprocessor_)
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
  CompressedSparsityPattern compressed_constrained_sparsity_pattern(n_dofs);
  DoFTools::make_sparsity_pattern(dof_handler,
    compressed_constrained_sparsity_pattern, constraints, false);
  constrained_sparsity_pattern.copy_from(compressed_constrained_sparsity_pattern);

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

  // set boundary indicators to distinguish incoming boundary
  setBoundaryIndicators();

  // determine Dirichlet nodes
  getDirichletNodes();
}

/**
 * Destructor.
 */
template<int dim>
Executioner<dim>::~Executioner()
{
}

/**
 * Assembles the inviscid steady-state matrix.
 */
template<int dim>
void Executioner<dim>::assembleInviscidSteadyStateMatrix()
{
  inviscid_ss_matrix = 0;

  // FE values, for assembly terms
  FEValues < dim
    > fe_values(fe, cell_quadrature,
      update_values | update_gradients | update_quadrature_points
        | update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // total cross section values at each quadrature point on cell
  std::vector<double> total_cross_section_values(n_q_points_cell);

  // cell iterator
  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active(), endc = dof_handler.end();
  // loop over cells
  unsigned int i_cell = 0;
  for (cell = dof_handler.begin_active(); cell != endc; ++cell, ++i_cell)
  {
    // initialize local matrix and rhs to zero
    cell_matrix = 0;

    // reinitialize FE values
    fe_values.reinit(cell);

    // get quadrature points on cell
    std::vector<Point<dim> > points(n_q_points_cell);
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
          cell_matrix(i, j) += (
          // divergence term
          fe_values[flux].value(i, q) * transport_direction
            * fe_values[flux].gradient(j, q) +
          // total interaction term
            fe_values[flux].value(i, q) * total_cross_section_values[q]
              * fe_values[flux].value(j, q)) * fe_values.JxW(q);
        } // end j
      } // end i
    } // end q

    // aggregate local matrix and rhs to global matrix and rhs
    constraints.distribute_local_to_global(cell_matrix, local_dof_indices,
      inviscid_ss_matrix);
  } // end cell
}

/**
 * Assembles the steady-state rhs.
 *
 * @param[in] t time at which to evaluate rhs
 */
template<int dim>
void Executioner<dim>::assembleSteadyStateRHS(const double &t)
{
  // reset steady-state rhs
  ss_rhs = 0;

  // set the time to be used in the source function
  source_function->set_time(t);

  // FE values, for assembly terms
  FEValues < dim
    > fe_values(fe, cell_quadrature,
      update_values | update_gradients | update_quadrature_points
        | update_JxW_values);

  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

  // source values at each quadrature point on cell
  std::vector<Point<dim> > points(n_q_points_cell);
  std::vector<double> source_values(n_q_points_cell);

  // cell iterator
  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active(), endc = dof_handler.end();
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
        cell_rhs(i) += fe_values[flux].value(i, q) * source_values[q]
          * fe_values.JxW(q);

    // aggregate local matrix and rhs to global matrix and rhs
    constraints.distribute_local_to_global(cell_rhs, local_dof_indices, ss_rhs);
  } // end cell
}

/**
 * Returns the final solution.
 */
template<int dim>
Vector<double> Executioner<dim>::getFinalSolution() const
{
  return new_solution;
}

/**
 * Sets the boundary indicators for each boundary face.
 *
 * The Dirichlet BC is applied only to the incoming boundary, so the transport
 * direction is compared against the normal vector of the face.
 */
template<int dim>
void Executioner<dim>::setBoundaryIndicators()
{
  // reset boundary indicators to zero
  typename DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active(), endc = dof_handler.end();
  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
    for (unsigned int face = 0; face < GeometryInfo < dim > ::faces_per_cell;
        ++face)
      if (cell->face(face)->at_boundary())
        cell->face(face)->set_boundary_indicator(0);

  // FE face values
  FEFaceValues<dim> fe_face_values(fe, face_quadrature,
    update_normal_vectors);
  // loop over cells
  for (cell = dof_handler.begin_active(); cell != endc; ++cell)
  {
    // loop over faces of cell
    for (unsigned int face = 0; face < GeometryInfo < dim > ::faces_per_cell;
        ++face)
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
          cell->face(face)->set_boundary_indicator(1);
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
template<int dim>
void Executioner<dim>::getDirichletNodes()
{
  // get map of Dirichlet dof indices to Dirichlet values
  std::map<unsigned int, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler, 1, ZeroFunction<dim>(),
      boundary_values);

  // extract dof indices from map
  dirichlet_nodes.clear();
  for (std::map<unsigned int, double>::iterator it = boundary_values.begin();
      it != boundary_values.end(); ++it)
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
template<int dim>
void Executioner<dim>::applyDirichletBC(
  SparseMatrix<double> & A,
  Vector<double>       & b,
  Vector<double>       & x,
  const double         & t)
{
   // set time for dirichlet function
   incoming_function->set_time(t);

   // create map of dofs to boundary values to be imposed
   std::map<unsigned int, double> boundary_values;
   VectorTools::interpolate_boundary_values(
     dof_handler,
     1,
     *incoming_function,
     boundary_values);

   // apply boundary values to system
   MatrixTools::apply_boundary_values(
     boundary_values,
     A,
     x,
     b);
}
