/**
 * \file Transport.cc
 * \brief Provides the function definitions for the Transport class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] params Transport equation parameters
 */
template <int dim>
Transport<dim>::Transport(const TransportParameters<dim> & params)
  : ConservationLaw<dim>(params),
    transport_parameters(params),
    extractor(0),
    cross_section_function(1),
    source_function(1)
{
}

template <int dim>
std::vector<std::string> Transport<dim>::get_component_names()
{
  std::vector<std::string> names(1, "angular_flux");
  return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Transport<
  dim>::get_component_interpretations()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      1, DataComponentInterpretation::component_is_scalar);

  return data_component_interpretation;
}

template <int dim>
void Transport<dim>::get_fe_extractors(
  std::vector<FEValuesExtractors::Scalar> & scalar_extractors,
  std::vector<FEValuesExtractors::Vector> &) const
{
  scalar_extractors.resize(1);
  scalar_extractors[0] = extractor;
}

template <int dim>
void Transport<dim>::define_problem()
{
  // determine problem name
  this->problem_name = transport_parameters.problem_name;

  // get path of source directory from #define created by CMake
  std::stringstream source_path_ss;
  source_path_ss << SOURCE_PATH;
  std::string source_path;
  source_path_ss >> source_path;

  // create problem parameters file name and determine if it exists
  std::string problem_parameters_file =
    source_path + "/problems/transport/" + this->problem_name;
  struct stat buffer;
  const bool file_exists = stat(problem_parameters_file.c_str(), &buffer) == 0;
  Assert(file_exists, ExcFileDoesNotExist(problem_parameters_file));

  // read problem parameters input file
  ParameterHandler parameter_handler;
  TransportProblemParameters<dim>::declare_parameters(parameter_handler);
  parameter_handler.read_input(problem_parameters_file);
  TransportProblemParameters<dim> problem_parameters;
  problem_parameters.get_parameters(parameter_handler);

  // assert number of dimensions is valid
  if (!problem_parameters.valid_in_1d)
    Assert(dim != 1, ExcImpossibleInDim(dim));
  if (!problem_parameters.valid_in_2d)
    Assert(dim != 2, ExcImpossibleInDim(dim));
  if (!problem_parameters.valid_in_3d)
    Assert(dim != 3, ExcImpossibleInDim(dim));

  // get unnormalized transport direction
  for (unsigned int d = 0; d < dim; ++d)
    if (d == 0)
      transport_direction[d] = problem_parameters.transport_direction_x;
    else if (d == 1)
      transport_direction[d] = problem_parameters.transport_direction_y;
    else
      transport_direction[d] = problem_parameters.transport_direction_z;

  // store constants for function parsers
  this->constants["incoming"] = problem_parameters.incoming_value;
  this->constants["sigma"] = problem_parameters.cross_section_value;
  this->constants["source"] = problem_parameters.source_value;
  this->constants["x_min"] = problem_parameters.x_start;

  // normalize transport direction
  double transport_direction_magnitude = 0.0;
  for (unsigned int d = 0; d < dim; ++d)
    transport_direction_magnitude += std::pow(transport_direction[d], 2);
  transport_direction_magnitude = std::sqrt(transport_direction_magnitude);
  for (unsigned int d = 0; d < dim; ++d)
    transport_direction[d] /= transport_direction_magnitude;

  // transport speed
  transport_speed = problem_parameters.transport_speed;

  // create string of variables to be used in function parser objects
  std::string variables = FunctionParser<dim>::default_variable_names();

  // add t for time-dependent function parser objects
  std::string time_dependent_variables = variables + ",t";

  // cross section
  cross_section_function.initialize(
    variables, problem_parameters.cross_section_string, this->constants, false);

  // source
  source_function.initialize(time_dependent_variables,
                             problem_parameters.source_string,
                             this->constants,
                             true);

  // domain
  if (problem_parameters.domain_shape == "hyper_cube")
  {
    const double x_start = problem_parameters.x_start;
    const double x_width = problem_parameters.x_width;
    this->domain_volume = std::pow(x_width, dim);
    GridGenerator::hyper_cube(this->triangulation, x_start, x_start + x_width);
  }
  else
  {
    Assert(false, ExcNotImplemented());
  }

  // set boundary indicators
  this->n_dirichlet_boundaries = 1;
  typename Triangulation<dim>::cell_iterator cell = this->triangulation.begin();
  typename Triangulation<dim>::cell_iterator endc = this->triangulation.end();
  for (; cell != endc; ++cell)
    for (unsigned int face = 0; face < this->faces_per_cell; ++face)
      if (cell->face(face)->at_boundary())
        cell->face(face)->set_boundary_id(0);

  // set boundary conditions type
  this->boundary_conditions_type = problem_parameters.boundary_conditions_type;

  // create boundary conditions
  if (this->boundary_conditions_type == "dirichlet")
  {
    auto derived_boundary_conditions =
      std::make_shared<DirichletBoundaryConditions<dim>>(this->fe,
                                                         this->face_quadrature);
    this->boundary_conditions = derived_boundary_conditions;

    // option to use exact solution as Dirichlet BC
    this->use_exact_solution_as_dirichlet_bc = false;

    // get Dirichlet function strings
    if (!this->use_exact_solution_as_dirichlet_bc)
    {
      this->dirichlet_function_strings.resize(this->n_dirichlet_boundaries);
      this->dirichlet_function_strings[0].resize(this->n_components);
      this->dirichlet_function_strings[0][0] =
        problem_parameters.dirichlet_function;
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }
  }
  else
  {
    Assert(false, ExcNotImplemented());
  }

  // initial conditions
  this->initial_conditions_strings[0] = problem_parameters.initial_condition;

  // exact solution
  if (problem_parameters.has_exact_solution)
  {
    this->has_exact_solution = true;

    if (problem_parameters.exact_solution_type == "function")
    {
      this->exact_solution_strings[0] = problem_parameters.exact_solution;

      // create and initialize function parser for exact solution
      std::shared_ptr<FunctionParser<dim>> exact_solution_function_derived =
        std::make_shared<FunctionParser<dim>>(this->parameters.n_components);
      exact_solution_function_derived->initialize(
        FunctionParser<dim>::default_variable_names() + ",t",
        this->exact_solution_strings,
        this->constants,
        true);
      this->exact_solution_function = exact_solution_function_derived;
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }
  }
  else
  {
    this->has_exact_solution = false;
  }

  // default end time
  this->has_default_end_time = problem_parameters.has_default_end_time;
  this->default_end_time = problem_parameters.default_end_time;
}

template <int dim>
void Transport<dim>::assemble_lumped_mass_matrix()
{
  FEValues<dim> fe_values(
    this->fe, this->cell_quadrature, update_values | update_JxW_values);

  std::vector<types::global_dof_index> local_dof_indices(this->dofs_per_cell);
  FullMatrix<double> local_mass(this->dofs_per_cell, this->dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    fe_values.reinit(cell);
    cell->get_dof_indices(local_dof_indices);

    local_mass = 0.0;

    // compute local contribution
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
        for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
        {
          local_mass(i, i) += fe_values[extractor].value(i, q) *
            fe_values[extractor].value(j, q) * fe_values.JxW(q);
        }

    // add to global mass matrix with contraints
    this->constraints.distribute_local_to_global(
      local_mass, local_dof_indices, this->lumped_mass_matrix_constrained);

    // add to global mass matrix without constraints
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
        this->lumped_mass_matrix.add(
          local_dof_indices[i], local_dof_indices[j], local_mass(i, j));
  }
}

/**
 * \brief Performs non-standard setup; calls function to distinguish boundary
 *        ID for incoming boundary.
 */
template <int dim>
void Transport<dim>::perform_nonstandard_setup()
{
  // distinguish boundary IDs for incoming boundary
  set_boundary_ids();
}

/**
 * \brief Sets the boundary IDs for each boundary face.
 *
 * The Dirichlet BC is applied only to the incoming boundary, so the transport
 * direction is compared against the normal vector of the face. The incoming
 * boundary is marked with a boundary ID of 0, while the rest of the boundary has
 * a boundary ID of 1.
 */
template <int dim>
void Transport<dim>::set_boundary_ids()
{
  // reset boundary indicators to zero
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
      if (cell->face(face)->at_boundary())
        cell->face(face)->set_boundary_id(1);

  // FE face values
  FEFaceValues<dim> fe_face_values(
    this->fe, this->face_quadrature, update_normal_vectors);

  // loop over cells
  for (cell = this->dof_handler.begin_active(); cell != endc; ++cell)
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
 * \brief Computes the steady-state flux vector.
 *
 * This function computes the steady-state flux vector \f$\mathrm{f}^n\f$:
 * \f[
 *   \mathrm{f}_i^n = \int\limits_{S_i}
 *     \varphi_i(\mathbf{x})\nabla\cdot\mathbf{f}(\tilde{\mathbf{u}}^n)dV \,.
 * \f]
 *
 * \param[in] dt time step size \f$\Delta t\f$
 * \param[in] solution solution vector \f$\mathrm{U}^n\f$
 * \param[out] ss_flux steady-state flux vector \f$\mathrm{f}^n\f$
 */
template <int dim>
void Transport<dim>::compute_ss_flux(const double & dt,
                                     const Vector<double> & solution,
                                     Vector<double> & ss_flux)
{
  // reset vector
  ss_flux = 0.0;

  FEValues<dim> fe_values(this->fe,
                          this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  Vector<double> cell_residual(this->dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);
  std::vector<double> solution_values(this->n_q_points_cell);
  std::vector<Tensor<1, dim>> solution_gradients(this->n_q_points_cell);

  // loop over cells
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reset cell residual
    cell_residual = 0;

    // reinitialize fe values for cell
    fe_values.reinit(cell);

    // get current solution values and gradients
    fe_values[extractor].get_function_values(solution, solution_values);
    fe_values[extractor].get_function_gradients(solution, solution_gradients);

    // get quadrature points on cell
    std::vector<Point<dim>> points(this->n_q_points_cell);
    points = fe_values.get_quadrature_points();

    // loop over quadrature points
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
      // compute cross section value at quadrature point
      const double cross_section_value = cross_section_function.value(points[q]);

      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        cell_residual(i) += transport_speed * fe_values[extractor].value(i, q) *
          (transport_direction * solution_gradients[q] +
           cross_section_value * solution_values[q]) *
          fe_values.JxW(q);
      }
    }

    // apply artificial diffusion
    /*
        this->artificial_diffusion->apply(
          this->viscosity, this->new_solution, cell, fe_values, cell_residual);
    */

    // apply boundary conditions
    this->boundary_conditions->apply(
      cell, fe_values, solution, dt, cell_residual);

    // aggregate local residual into global residual
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      ss_flux(local_dof_indices[i]) += cell_residual(i);
  }
}

/**
 * \brief Computes the steady-state right hand side vector.
 *
 * This function computes the steady-state flux vector \f$\mathrm{q}^n\f$:
 * \f[
 *   \mathrm{q}_i^n = v \int\limits_{S_i}
 *     \varphi_i(\mathbf{x}) q(t^n,\mathbf{x})dV \,,
 * \f]
 * where \f$v\f$ is the transport speed.
 *
 * \param[in] t time at which to evaluate residual \f$t^n\f$
 * \param[out] ss_flux steady-state rhs vector \f$\mathrm{q}^n\f$
 */
template <int dim>
void Transport<dim>::compute_ss_rhs(const double & t, Vector<double> & ss_rhs)
{
  // reset vector
  ss_rhs = 0.0;

  FEValues<dim> fe_values(this->fe,
                          this->cell_quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  Vector<double> cell_residual(this->dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);

  // set time for source function
  source_function.set_time(t);

  // loop over cells
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reset cell residual
    cell_residual = 0;

    // reinitialize fe values for cell
    fe_values.reinit(cell);

    // get quadrature points on cell
    std::vector<Point<dim>> points(this->n_q_points_cell);
    points = fe_values.get_quadrature_points();

    // loop over quadrature points
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
      // compute source value at quadrature point
      const double source_value = source_function.value(points[q]);

      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        cell_residual(i) += transport_speed * fe_values[extractor].value(i, q) *
          source_value * fe_values.JxW(q);
      }
    }

    // aggregate local residual into global residual
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      ss_rhs(local_dof_indices[i]) += cell_residual(i);
  }
}

template <int dim>
void Transport<dim>::update_flux_speeds()
{
  // loop over cells to update max flux speed for each cell
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
    this->max_flux_speed_cell[cell] = transport_speed;

  // reset max flux speed
  this->max_flux_speed = transport_speed;
}

/**
 * \brief Creates an entropy object.
 *
 * \return pointer to created entropy object
 */
template <int dim>
std::shared_ptr<Entropy<dim>> Transport<dim>::create_entropy() const
{
  auto entropy = std::make_shared<ScalarEntropy<dim>>(this->domain_volume,
                                                      this->dof_handler,
                                                      this->fe,
                                                      this->triangulation,
                                                      this->cell_quadrature,
                                                      this->face_quadrature);
  return entropy;
}

/**
 * \brief Creates a max wave speed object
 *
 * \return pointer to created max wave speed object
 */
template <int dim>
std::shared_ptr<MaxWaveSpeed<dim>> Transport<dim>::create_max_wave_speed() const
{
  auto max_wave_speed =
    std::make_shared<ConstantMaxWaveSpeed<dim>>(transport_speed);
  return max_wave_speed;
}
