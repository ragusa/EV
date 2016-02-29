/**
 * \file Burgers.cc
 * \brief Provides the function definitions for the Burgers class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] params Burgers equation parameters
 */
template <int dim>
Burgers<dim>::Burgers(const BurgersParameters<dim> & params)
  : ConservationLaw<dim>(params),
    burgers_parameters(params),
    velocity_extractor(0)
{
}

template <int dim>
std::vector<std::string> Burgers<dim>::get_component_names()
{
  std::vector<std::string> names(1, "velocity");
  return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Burgers<
  dim>::get_component_interpretations()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      1, DataComponentInterpretation::component_is_scalar);

  return data_component_interpretation;
}

template <int dim>
void Burgers<dim>::get_fe_extractors(
  std::vector<FEValuesExtractors::Scalar> & scalar_extractors,
  std::vector<FEValuesExtractors::Vector> &) const
{
  scalar_extractors.resize(1);
  scalar_extractors[0] = velocity_extractor;
}

template <int dim>
void Burgers<dim>::define_problem()
{
  // determine problem name
  this->problem_name = burgers_parameters.problem_name;

  // get path of source directory from #define created by CMake
  std::stringstream source_path_ss;
  source_path_ss << SOURCE_PATH;
  std::string source_path;
  source_path_ss >> source_path;

  // create problem parameters file name and determine if it exists
  std::string problem_parameters_file =
    source_path + "/problems/burgers/" + this->problem_name;
  struct stat buffer;
  const bool file_exists = stat(problem_parameters_file.c_str(), &buffer) == 0;
  Assert(file_exists, ExcFileDoesNotExist(problem_parameters_file));

  // read problem parameters input file
  ParameterHandler parameter_handler;
  BurgersProblemParameters<dim>::declare_parameters(parameter_handler);
  parameter_handler.read_input(problem_parameters_file);
  BurgersProblemParameters<dim> problem_parameters;
  problem_parameters.get_parameters(parameter_handler);

  // assert number of dimensions is valid
  if (!problem_parameters.valid_in_1d)
    Assert(dim != 1, ExcImpossibleInDim(dim));
  if (!problem_parameters.valid_in_2d)
    Assert(dim != 2, ExcImpossibleInDim(dim));
  if (!problem_parameters.valid_in_3d)
    Assert(dim != 3, ExcImpossibleInDim(dim));

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

  // constants
  const double x_interface = problem_parameters.x_interface;
  const double u_left = problem_parameters.u_left;
  const double u_right = problem_parameters.u_right;
  this->constants["x_interface"] = x_interface;
  this->constants["u_left"] = u_left;
  this->constants["u_right"] = u_right;

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
  this->initial_conditions_strings[0] = problem_parameters.initial_conditions;

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
void Burgers<dim>::assemble_lumped_mass_matrix()
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
          local_mass(i, i) += fe_values[velocity_extractor].value(i, q) *
            fe_values[velocity_extractor].value(j, q) * fe_values.JxW(q);
        }

    // add to global mass matrix with contraints
    this->constraints.distribute_local_to_global(
      local_mass, local_dof_indices, this->lumped_mass_matrix);
  }
}

/**
 * \brief Computes the steady-state residual.
 *
 * Rearranging the Burgers equation,
 * \f[
 *   \frac{\partial u}{\partial t}
 *   = - u\mathbf{v}\cdot\nabla u .
 * \f]
 * Substituting the approximate FEM solution and testing with a test function
 * \f$\varphi_i\f$ gives the weak form for degree of freedom \f$i\f$:
 * \f[
 *   \left(\varphi_i,\frac{\partial u_h}{\partial t}\right)_\Omega
 *   = - \left(\varphi_i,u_h\mathbf{v}\cdot\nabla u_h\right)_\Omega .
 * \f]
 * This yields a discrete system
 * \f[
 *   \mathbf{M}\frac{d\mathbf{U}}{dt} = \mathbf{r} ,
 * \f]
 * where \f$\mathbf{M}\f$ is the mass matrix and the steady-state residual
 * \f$\mathbf{r}\f$ is given by
 * \f[
 *   r_i = - \left(\varphi_i,u_h\mathbf{v}\cdot\nabla u_h\right)_\Omega
 * \f]
 *
 *  \param[in] t time at which to evaluate residual \f$t\f$
 *  \param[in] dt time step size \f$\Delta t\f$
 *  \param[out] r steady-state residual \f$\mathbf{r}\f$
 */
template <int dim>
void Burgers<dim>::compute_ss_residual(const double &,
                                       const double & dt,
                                       Vector<double> & f)
{
  // reset vector
  f = 0.0;

  FEValues<dim> fe_values(this->fe,
                          this->cell_quadrature,
                          update_values | update_gradients | update_JxW_values);

  Vector<double> cell_residual(this->dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);
  std::vector<double> solution_values(this->n_q_points_cell);
  std::vector<Tensor<1, dim>> solution_gradients(this->n_q_points_cell);
  std::vector<Tensor<1, dim>> dfdu(this->n_q_points_cell);

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
    fe_values[velocity_extractor].get_function_values(this->new_solution,
                                                      solution_values);
    fe_values[velocity_extractor].get_function_gradients(this->new_solution,
                                                         solution_gradients);

    // loop over quadrature points
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
      // compute derivative of flux
      for (int d = 0; d < dim; ++d)
        dfdu[q][d] = solution_values[q];

      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        cell_residual(i) += -fe_values[velocity_extractor].value(i, q) * dfdu[q] *
          solution_gradients[q] * fe_values.JxW(q);
      }
    }

    // apply artificial diffusion
    this->artificial_diffusion->apply(
      this->viscosity, this->new_solution, cell, fe_values, cell_residual);

    // apply boundary conditions
    this->boundary_conditions->apply(
      cell, fe_values, this->new_solution, dt, cell_residual);

    // aggregate local residual into global residual
    cell->get_dof_indices(local_dof_indices);
    this->constraints.distribute_local_to_global(
      cell_residual, local_dof_indices, f);
  } // end cell loop

  // if artificial diffusion is algebraic, then apply it here
  this->artificial_diffusion->apply_algebraic_diffusion(this->new_solution, f);
}

template <int dim>
void Burgers<dim>::update_flux_speeds()
{
  FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);
  Tensor<1, dim> dfdu;
  std::vector<double> velocity(this->n_q_points_cell);

  // reset max flux speed
  this->max_flux_speed = 0.0;

  // loop over cells to compute first order viscosity at each quadrature point
  typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                          .begin_active(),
                                                 endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    fe_values.reinit(cell);
    fe_values[velocity_extractor].get_function_values(this->new_solution,
                                                      velocity);

    this->max_flux_speed_cell[cell] = 0.0;
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
      for (unsigned int d = 0; d < dim; ++d)
        dfdu[d] = velocity[q];
      this->max_flux_speed_cell[cell] =
        std::max(this->max_flux_speed_cell[cell], dfdu.norm());
    }

    // get max flux speed
    this->max_flux_speed =
      std::max(this->max_flux_speed, this->max_flux_speed_cell[cell]);
  }
}

/**
 * \brief Creates an entropy object.
 *
 * \return pointer to created entropy object
 */
template <int dim>
std::shared_ptr<Entropy<dim>> Burgers<dim>::create_entropy() const
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
std::shared_ptr<MaxWaveSpeed<dim>> Burgers<dim>::create_max_wave_speed() const
{
  auto max_wave_speed = std::make_shared<BurgersMaxWaveSpeed<dim>>();
  return max_wave_speed;
}
