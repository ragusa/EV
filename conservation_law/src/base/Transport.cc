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
Transport<dim>::Transport(const TransportRunParameters & params)
  : ConservationLaw<dim>(params, 1, true),
    transport_parameters(params),
    problem_parameters(params.problem_name, false),
    extractor(0)
{
  // point base class problem parameters to derived class problem parameters
  this->problem_base_parameters = &problem_parameters;
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
  // get path of source directory from #define created by CMake
  std::stringstream source_path_ss;
  source_path_ss << SOURCE_PATH;
  std::string source_path;
  source_path_ss >> source_path;

  // create problem parameters file name and determine if it exists
  std::string problem_parameters_file =
    source_path + "/problems/transport/" + transport_parameters.problem_name;

  // get and process the problem parameters
  problem_parameters.get_and_process_parameters(problem_parameters_file,
                                                this->triangulation,
                                                this->fe,
                                                this->face_quadrature);
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

    // add to global mass matrix
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
}

/**
 * \brief Computes the inviscid steady-state matrix \f$\mathbf{A}\f$.
 *
 * This function computes the inviscid steady-state matrix \f$\mathbf{A}\f$:
 * \f[
 *   A_{i,j} = \int\limits_{S_{i,j}}
 *     \left(\mathbf{f}'(\tilde{\mathbf{u}}^n)\cdot\nabla\varphi_j(\mathbf{x})
 *       + \sigma(\mathbf{x})\varphi_j(\mathbf{x})\right)\varphi_i(\mathbf{x})
 *     dV \,.
 * \f]
 *
 * \param[in] solution  solution vector \f$\mathrm{U}^n\f$
 * \param[out] matrix   inviscid steady-state matrix \f$\mathrm{A}\f$
 */
template <int dim>
void Transport<dim>::compute_inviscid_ss_matrix(const Vector<double> & solution,
                                                SparseMatrix<double> & matrix)
{
  // reset
  matrix = 0.0;

  FEValues<dim> fe_values(this->fe,
                          this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_matrix(this->dofs_per_cell, this->dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);
  std::vector<Tensor<1, dim>> solution_gradients(this->n_q_points_cell);

  // loop over cells
  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reset cell matrix
    cell_matrix = 0;

    // reinitialize fe values for cell
    fe_values.reinit(cell);

    // get quadrature points on cell
    std::vector<Point<dim>> points(this->n_q_points_cell);
    points = fe_values.get_quadrature_points();

    // loop over quadrature points
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
      // compute cross section value at quadrature point
      const double cross_section_value =
        problem_parameters.cross_section_function.value(points[q]);

      // loop over test function i
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        // loop over test function j
        for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
        {
          cell_matrix(i, j) +=
            (problem_parameters.transport_direction *
               fe_values[extractor].gradient(j, q) +
             cross_section_value * fe_values[extractor].value(j, q)) *
            problem_parameters.transport_speed *
            fe_values[extractor].value(i, q) * fe_values.JxW(q);
        }
      }
    }

    // apply boundary conditions
    problem_parameters.boundary_conditions->apply(
      cell, fe_values, solution, 0.0, cell_matrix);

    // aggregate local residual into global residual
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
        matrix.add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
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
      const double cross_section_value =
        problem_parameters.cross_section_function.value(points[q]);

      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        cell_residual(i) += problem_parameters.transport_speed *
          fe_values[extractor].value(i, q) *
          (problem_parameters.transport_direction * solution_gradients[q] +
           cross_section_value * solution_values[q]) *
          fe_values.JxW(q);
      }
    }

    // apply boundary conditions
    problem_parameters.boundary_conditions->apply(
      cell, fe_values, solution, dt, cell_residual);

    // aggregate local residual into global residual
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      ss_flux(local_dof_indices[i]) += cell_residual(i);
  }
}

/**
 * \brief Computes the steady-state reaction vector.
 *
 * This function computes the steady-state reaction vector, which has the
 * entries
 * \f[
 *   \sigma_i = \int\limits_{S_i}\varphi_i(\mathbf{x})\sigma(\mathbf{x})dV
 *     = \sum\limits_j\int\limits_{S_{i,j}}
 *     \varphi_i(\mathbf{x})\varphi_j(\mathbf{x})\sigma(\mathbf{x})dV \,,
 * \f]
 * where \f$\sigma(\mathbf{x})\f$ is the reaction coefficient of the
 * conservation law equation when it is put in the form
 * \f[
 *   \frac{\partial\mathbf{u}}{\partial t}
 *   + \nabla \cdot \mathbf{F}(\mathbf{u}) + \sigma(\mathbf{x})\mathbf{u}
 *   = \mathbf{0} \,.
 * \f]
 *
 * \param[out] ss_reaction steady-state reaction vector
 */
template <int dim>
void Transport<dim>::compute_ss_reaction(Vector<double> & ss_reaction)
{
  // reset vector
  ss_reaction = 0.0;

  FEValues<dim> fe_values(this->fe,
                          this->cell_quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);

  Vector<double> cell_residual(this->dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);

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
      // compute cross section value at quadrature point
      const double cross_section_value =
        problem_parameters.cross_section_function.value(points[q]);

      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        cell_residual(i) += problem_parameters.transport_speed *
          fe_values[extractor].value(i, q) * cross_section_value *
          fe_values.JxW(q);
      }
    }

    // aggregate local residual into global residual
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      ss_reaction(local_dof_indices[i]) += cell_residual(i);
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
  problem_parameters.source_function.set_time(t);

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
      const double source_value =
        problem_parameters.source_function.value(points[q]);

      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        cell_residual(i) += problem_parameters.transport_speed *
          fe_values[extractor].value(i, q) * source_value * fe_values.JxW(q);
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
    this->max_flux_speed_cell[cell] = problem_parameters.transport_speed;

  // reset max flux speed
  this->max_flux_speed = problem_parameters.transport_speed;
}

/**
 * \brief Creates an entropy object.
 *
 * \return pointer to created entropy object
 */
template <int dim>
std::shared_ptr<Entropy<dim>> Transport<dim>::create_entropy() const
{
  auto entropy =
    std::make_shared<TransportEntropy<dim>>(problem_parameters.domain_volume,
                                            this->dof_handler,
                                            this->fe,
                                            this->cell_quadrature,
                                            this->face_quadrature,
                                            problem_parameters);
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
  auto max_wave_speed = std::make_shared<ConstantMaxWaveSpeed<dim>>(
    problem_parameters.transport_speed);
  return max_wave_speed;
}

/**
 * \brief Creates an FCT object and returns the pointer.
 *
 * \return pointer to new FCT object
 */
template <int dim>
std::shared_ptr<ExplicitEulerFCT<dim>> Transport<dim>::create_fct() const
{
  // determine if star states are to be used in FCT bounds
  const bool use_star_states_in_fct_bounds =
    this->are_star_states && this->parameters.use_star_states_in_fct_bounds;
  Assert(!use_star_states_in_fct_bounds, ExcNotImplemented());

  // create FCT object
  auto fct =
    std::make_shared<TransportExplicitEulerFCT<dim>>(transport_parameters,
                                                     problem_parameters,
                                                     this->dof_handler,
                                                     this->fe,
                                                     this->dirichlet_dof_indices,
                                                     this->consistent_mass_matrix,
                                                     this->lumped_mass_matrix);

  return fct;
}
