/**
 * \file ShallowWater.cc
 * \brief Provides the function definitions for the ShallowWater class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] params shallow water equation parameters
 */
template <int dim>
ShallowWater<dim>::ShallowWater(const ShallowWaterParameters<dim> & params)
  : ConservationLaw<dim>(params, dim + 1, true),
    sw_parameters(params),
    height_extractor(0),
    momentum_extractor(1),
    fe_bathymetry(params.degree),
    dof_handler_bathymetry(this->triangulation)
{
  // shallow water equations cannot be 3-D
  Assert(dim < 3, ExcImpossibleInDim(dim));
}

template <int dim>
std::vector<std::string> ShallowWater<dim>::get_component_names()
{
  std::vector<std::string> names(1 + dim);

  names[0] = "height";
  for (int d = 0; d < dim; ++d)
    names[1 + d] = "momentum";

  return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
  ShallowWater<dim>::get_component_interpretations()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretations(dim + 1);

  component_interpretations[0] = DataComponentInterpretation::component_is_scalar;
  for (int d = 0; d < dim; ++d)
    component_interpretations[1 + d] =
      DataComponentInterpretation::component_is_part_of_vector;

  return component_interpretations;
}

template <int dim>
void ShallowWater<dim>::get_fe_extractors(
  std::vector<FEValuesExtractors::Scalar> & scalar_extractors,
  std::vector<FEValuesExtractors::Vector> & vector_extractors) const
{
  scalar_extractors.resize(1);
  scalar_extractors[0] = height_extractor;

  vector_extractors.resize(1);
  vector_extractors[0] = momentum_extractor;
}

template <int dim>
void ShallowWater<dim>::define_problem()
{
  // determine problem name
  this->problem_name = sw_parameters.problem_name;

  // get path of source directory from #define created by CMake
  std::stringstream source_path_ss;
  source_path_ss << SOURCE_PATH;
  std::string source_path;
  source_path_ss >> source_path;

  // create problem parameters file name and determine if it exists
  std::string problem_parameters_file =
    source_path + "/problems/shallowwater/" + this->problem_name;
  struct stat buffer;
  const bool file_exists = stat(problem_parameters_file.c_str(), &buffer) == 0;
  Assert(file_exists, ExcFileDoesNotExist(problem_parameters_file));

  // read problem parameters input file
  ParameterHandler parameter_handler;
  ShallowWaterProblemParameters<dim>::declare_parameters(parameter_handler);
  parameter_handler.read_input(problem_parameters_file);
  ShallowWaterProblemParameters<dim> problem_parameters;
  problem_parameters.get_parameters(parameter_handler);

  // assert number of dimensions is valid
  if (!problem_parameters.valid_in_1d)
    Assert(dim != 1, ExcImpossibleInDim(dim));
  if (!problem_parameters.valid_in_2d)
    Assert(dim != 2, ExcImpossibleInDim(dim));

  // domain
  if (problem_parameters.domain_shape == "hyper_cube")
  {
    const double x_start = problem_parameters.x_start;
    const double x_width = problem_parameters.x_width;
    this->domain_volume = std::pow(x_width, dim);
    GridGenerator::hyper_cube(this->triangulation, x_start, x_start + x_width);
  }
  else if (problem_parameters.domain_shape == "hyper_rectangle")
  {
    Point<dim> point_start;
    Point<dim> point_end;
    const double x_start = problem_parameters.x_start;
    const double x_width = problem_parameters.x_width;
    point_start[0] = x_start;
    point_end[0] = x_start + x_width;
    if (dim == 1)
    {
      this->domain_volume = x_width;
    }
    else
    {
      const double y_start = problem_parameters.y_start;
      const double y_width = problem_parameters.y_width;
      point_start[1] = y_start;
      point_end[1] = y_start + y_width;
      this->domain_volume = x_width * y_width;
    }
    GridGenerator::hyper_rectangle(this->triangulation, point_start, point_end);
  }
  else
  {
    Assert(false, ExcNotImplemented());
  }

  // constants
  gravity = problem_parameters.gravity;
  const double x_interface = problem_parameters.x_interface;
  const double h_left = problem_parameters.h_left;
  const double h_right = problem_parameters.h_right;
  const double h_unperturbed = problem_parameters.h_unperturbed;
  const double h_perturbed = problem_parameters.h_perturbed;
  const double u_left = problem_parameters.u_left;
  const double u_right = problem_parameters.u_right;
  const double bump_height = problem_parameters.bump_height;
  const double bump_x_center = problem_parameters.bump_x_center;
  const double bump_y_center = problem_parameters.bump_y_center;
  const double bump_x_width = problem_parameters.bump_x_width;
  const double bump_y_width = problem_parameters.bump_y_width;
  const double perturbation_x_center = problem_parameters.perturbation_x_center;
  const double perturbation_y_center = problem_parameters.perturbation_y_center;
  const double perturbation_x_width = problem_parameters.perturbation_x_width;
  const double perturbation_y_width = problem_parameters.perturbation_y_width;
  this->constants["x_interface"] = x_interface;
  this->constants["h_left"] = h_left;
  this->constants["h_right"] = h_right;
  this->constants["h_unperturbed"] = h_unperturbed;
  this->constants["h_perturbed"] = h_perturbed;
  this->constants["u_left"] = u_left;
  this->constants["u_right"] = u_right;
  this->constants["bump_height"] = bump_height;
  this->constants["bump_x_center"] = bump_x_center;
  this->constants["bump_y_center"] = bump_y_center;
  this->constants["bump_x_width"] = bump_x_width;
  this->constants["bump_y_width"] = bump_y_width;
  this->constants["bump_left"] = bump_x_center - 0.5 * bump_x_width;
  this->constants["bump_right"] = bump_x_center + 0.5 * bump_x_width;
  this->constants["perturbation_x_center"] = perturbation_x_center;
  this->constants["perturbation_y_center"] = perturbation_y_center;
  this->constants["perturbation_x_width"] = perturbation_x_width;
  this->constants["perturbation_y_width"] = perturbation_y_width;

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
        problem_parameters.dirichlet_function_height;
      this->dirichlet_function_strings[0][1] =
        problem_parameters.dirichlet_function_momentumx;
      if (dim == 2)
        this->dirichlet_function_strings[0][2] =
          problem_parameters.dirichlet_function_momentumy;
    }
    else
    {
      Assert(false, ExcNotImplemented());
    }
  }
  else if (this->boundary_conditions_type == "none")
  {
    std::shared_ptr<ShallowWaterNoBC<dim>> derived_boundary_conditions =
      std::make_shared<ShallowWaterNoBC<dim>>(
        this->fe, this->face_quadrature, gravity);
    this->boundary_conditions = derived_boundary_conditions;
  }
  else if (this->boundary_conditions_type == "characteristic_open")
  {
    std::shared_ptr<ShallowWaterSubcriticalOpenBC1D<dim>>
      derived_boundary_conditions =
        std::make_shared<ShallowWaterSubcriticalOpenBC1D<dim>>(
          this->fe, this->face_quadrature, gravity, h_unperturbed, h_unperturbed);
    this->boundary_conditions = derived_boundary_conditions;
  }
  else if (this->boundary_conditions_type == "wall")
  {
    std::shared_ptr<ShallowWaterWallBC<dim>> derived_boundary_conditions =
      std::make_shared<ShallowWaterWallBC<dim>>(
        this->fe, this->face_quadrature, gravity);
    this->boundary_conditions = derived_boundary_conditions;
  }
  else if (this->boundary_conditions_type == "characteristic_wall")
  {
    std::shared_ptr<ShallowWaterSubcriticalWallBC1D<dim>>
      derived_boundary_conditions =
        std::make_shared<ShallowWaterSubcriticalWallBC1D<dim>>(
          this->fe, this->face_quadrature, gravity);
    this->boundary_conditions = derived_boundary_conditions;
  }
  else
  {
    Assert(false, ExcNotImplemented());
  }

  // initial conditions
  this->initial_conditions_strings[0] =
    problem_parameters.initial_conditions_height;
  this->initial_conditions_strings[1] =
    problem_parameters.initial_conditions_momentumx;
  if (dim == 2)
    this->initial_conditions_strings[2] =
      problem_parameters.initial_conditions_momentumy;

  // exact solution
  if (problem_parameters.has_exact_solution)
  {
    this->has_exact_solution = true;

    if (problem_parameters.exact_solution_type == "riemann")
    {
      // create and initialize Riemann solver for exact solution
      std::shared_ptr<ShallowWaterRiemannSolver<dim>>
        exact_solution_function_derived =
          std::make_shared<ShallowWaterRiemannSolver<dim>>(
            h_left, u_left, h_right, u_right, gravity, x_interface);
      this->exact_solution_function = exact_solution_function_derived;
    }
    else if (problem_parameters.exact_solution_type == "function")
    {
      this->exact_solution_strings[0] = problem_parameters.exact_solution_height;
      this->exact_solution_strings[1] =
        problem_parameters.exact_solution_momentumx;
      if (dim == 2)
        this->exact_solution_strings[2] =
          problem_parameters.exact_solution_momentumy;

      // create and initialize function parser for exact solution
      std::shared_ptr<FunctionParser<dim>> exact_solution_function_derived =
        std::make_shared<FunctionParser<dim>>(this->n_components);
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

  // initialize bathymetry function
  std::shared_ptr<FunctionParser<dim>> bathymetry_function_derived =
    std::make_shared<FunctionParser<dim>>();
  std::map<std::string, double> constants;
  bathymetry_function_derived->initialize(
    FunctionParser<dim>::default_variable_names(),
    problem_parameters.bathymetry_function,
    this->constants,
    false);
  bathymetry_function = bathymetry_function_derived;

  // default end time
  this->has_default_end_time = problem_parameters.has_default_end_time;
  this->default_end_time = problem_parameters.default_end_time;
}

/**
 * \brief Creates an entropy object.
 *
 * \return pointer to created entropy object
 */
template <int dim>
std::shared_ptr<Entropy<dim>> ShallowWater<dim>::create_entropy() const
{
  auto entropy =
    std::make_shared<ShallowWaterEntropy<dim>>(sw_parameters,
                                               height_extractor,
                                               momentum_extractor,
                                               gravity,
                                               bathymetry_vector,
                                               this->domain_volume,
                                               this->dof_handler,
                                               this->fe,
                                               this->triangulation,
                                               this->cell_quadrature,
                                               this->face_quadrature);
  return entropy;
}

/**
 * \brief Creates a star state object.
 *
 * \return pointer to created star state object
 */
template <int dim>
std::shared_ptr<StarState<dim>> ShallowWater<dim>::create_star_state() const
{
  auto star_state =
    std::make_shared<ShallowWaterStarState<dim>>(gravity,
                                                 this->gradient_matrix,
                                                 this->dof_handler,
                                                 this->triangulation,
                                                 this->n_components);

  return star_state;
}

/**
 * \brief Creates a max wave speed object
 *
 * \return pointer to created max wave speed object
 */
template <int dim>
std::shared_ptr<MaxWaveSpeed<dim>> ShallowWater<dim>::create_max_wave_speed()
  const
{
  auto max_wave_speed =
    std::make_shared<ShallowWaterMaxWaveSpeed<dim>>(this->star_state, gravity);

  return max_wave_speed;
}

/**
 * \brief Interpolates the bathymetry FE vector from its function.
 */
template <int dim>
void ShallowWater<dim>::perform_nonstandard_setup()
{
  // interpolate bathymetry vector
  dof_handler_bathymetry.clear();
  dof_handler_bathymetry.distribute_dofs(fe_bathymetry);
  bathymetry_vector.reinit(dof_handler_bathymetry.n_dofs());
  VectorTools::interpolate(
    dof_handler_bathymetry, *bathymetry_function, bathymetry_vector);
}

template <int dim>
void ShallowWater<dim>::assemble_lumped_mass_matrix()
{
  FEValues<dim> fe_values(
    this->fe, this->cell_quadrature, update_values | update_JxW_values);

  std::vector<types::global_dof_index> local_dof_indices(this->dofs_per_cell);
  FullMatrix<double> local_mass(this->dofs_per_cell, this->dofs_per_cell);

  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
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
          local_mass(i, i) += (fe_values[height_extractor].value(i, q) *
                                 fe_values[height_extractor].value(j, q) +
                               fe_values[momentum_extractor].value(i, q) *
                                 fe_values[momentum_extractor].value(j, q)) *
            fe_values.JxW(q);
        }

    // add to global mass matrix
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
        this->lumped_mass_matrix.add(
          local_dof_indices[i], local_dof_indices[j], local_mass(i, j));
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
 * Rearranging the continuity equation, substituting the approximate FEM
 * solution and testing with a test function \f$\varphi_i^h\f$ gives its
 * weak form for degree of freedom \f$i\f$:
 * \f[
 *   \left(\varphi_i^h,\frac{\partial h_h}{\partial t}\right)_\Omega
 *   = - \left(\varphi_i^h,\nabla\cdot\mathbf{q}_h\right)_\Omega .
 * \f]
 * Integrating by parts gives
 * \f[
 *   \left(\varphi_i^h,\frac{\partial h_h}{\partial t}\right)_\Omega
 *   = \left(\nabla\varphi_i^h,\mathbf{q}_h\right)_\Omega
 *   - \left(\varphi_i^h,\mathbf{q}_h\cdot\mathbf{n}\right)_{\partial\Omega} .
 * \f]
 * Rearranging the momentum equation, substituting the approximate FEM
 * solution and testing with a test function \f$\varphi_i^{\mathbf{q}}\f$
 * gives its weak form for degree of freedom \f$i\f$:
 * \f[
 *   \left(\varphi_i^{\mathbf{q}},\frac{\partial\mathbf{q}_h}{\partial t}
 *   \right)_\Omega
 *   = - \left(\varphi_i^{\mathbf{q}},\nabla\cdot(\mathbf{q}\otimes\mathbf{v}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\right)_\Omega
 *   - \left(\varphi_i^{\mathbf{q}},g h_h\nabla b\right)_\Omega .
 * \f]
 * Integrating by parts gives
 * \f[
 *   \left(\varphi_i^{\mathbf{q}},\frac{\partial\mathbf{q}_h}{\partial t}
 *   \right)_\Omega
 *   = \left(\nabla\varphi_i^{\mathbf{q}},(\mathbf{q}\otimes\mathbf{v}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\right)_\Omega
 *   - \left(\varphi_i^{\mathbf{q}},\mathbf{q}\otimes\mathbf{v}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\cdot\mathbf{n}\right)_{\partial\Omega}
 *   - \left(\varphi_i^{\mathbf{q}},g h_h\nabla b\right)_\Omega .
 * \f]
 * This yields a discrete system
 * \f[
 *   \mathbf{M}\frac{d\mathbf{U}}{dt} = \mathbf{r} ,
 * \f]
 * where \f$\mathbf{M}\f$ is the mass matrix and the steady-state residual
 * \f$\mathbf{r}\f$ is given by
 * \f[
 *   r_i =
 *   \left(\nabla\varphi_i^h,\mathbf{q}_h\right)_\Omega
 *   - \left(\varphi_i^h,\mathbf{q}_h\cdot\mathbf{n}\right)_{\partial\Omega}
 *   + \left(\nabla\varphi_i^{\mathbf{q}},(\mathbf{q}\otimes\mathbf{v}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\right)_\Omega
 *   - \left(\varphi_i^{\mathbf{q}},(\mathbf{q}\otimes\mathbf{v}
 *   + \frac{1}{2}g h^2\mathbf{I})_h\cdot\mathbf{n}\right)_{\partial\Omega}
 *   - \left(\varphi_i^{\mathbf{q}},g h_h\nabla b\right)_\Omega .
 * \f]
 *
 * \param[in] dt time step size \f$\Delta t\f$
 * \param[in] solution solution vector \f$\mathrm{U}^n\f$
 * \param[out] ss_flux steady-state flux vector \f$\mathrm{f}^n\f$
 */
template <int dim>
void ShallowWater<dim>::compute_ss_flux(const double & dt,
                                        const Vector<double> & solution,
                                        Vector<double> & ss_flux)
{
  // reset vector
  ss_flux = 0;

  FEValues<dim> fe_values(this->fe,
                          this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FEValues<dim> fe_values_bathymetry(
    fe_bathymetry, this->cell_quadrature, update_gradients);

  Vector<double> cell_residual(this->dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);

  // loop over cells
  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end(),
       cell_bathymetry = dof_handler_bathymetry.begin_active();
  for (; cell != endc; ++cell, ++cell_bathymetry)
  {
    // reset cell residual
    cell_residual = 0;

    // reinitialize fe values for cell
    fe_values.reinit(cell);
    fe_values_bathymetry.reinit(cell_bathymetry);

    // get current solution values
    std::vector<double> height(this->n_q_points_cell);
    std::vector<Tensor<1, dim>> momentum(this->n_q_points_cell);
    fe_values[height_extractor].get_function_values(solution, height);
    fe_values[momentum_extractor].get_function_values(solution, momentum);

    // get solution gradients
    std::vector<Tensor<1, dim>> height_gradient(this->n_q_points_cell);
    std::vector<Tensor<2, dim>> momentum_gradient(this->n_q_points_cell);
    fe_values[height_extractor].get_function_gradients(solution, height_gradient);
    fe_values[momentum_extractor].get_function_gradients(solution,
                                                         momentum_gradient);

    // compute gradients of bathymetry function
    std::vector<Tensor<1, dim>> bathymetry_gradient(this->n_q_points_cell);
    fe_values_bathymetry.get_function_gradients(bathymetry_vector,
                                                bathymetry_gradient);

    // compute inviscid fluxes
    std::vector<Tensor<1, dim>> height_inviscid_flux(this->n_q_points_cell);
    std::vector<Tensor<2, dim>> momentum_inviscid_flux(this->n_q_points_cell);
    compute_inviscid_fluxes(
      height, momentum, height_inviscid_flux, momentum_inviscid_flux);

    // loop over quadrature points
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        cell_residual(i) +=
          (
            // height flux
            -fe_values[height_extractor].gradient(i, q) * height_inviscid_flux[q]
            // momentum flux
            -
            double_contract<0, 0, 1, 1>(
              fe_values[momentum_extractor].gradient(i, q),
              momentum_inviscid_flux[q])
            // bathymetry source term
            +
            fe_values[momentum_extractor].value(i, q) * gravity * height[q] *
              bathymetry_gradient[q]) *
          fe_values.JxW(q);
      }
    }

    // apply artificial diffusion
    /*
        this->artificial_diffusion->apply(
          this->viscosity, solution, cell, fe_values, cell_residual);
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
 *   \mathrm{q}_i^n = \int\limits_{S_i}
 *     \varphi_i(\mathbf{x}) q(t^n,\mathbf{x})dV \,.
 * \f]
 *
 * \param[in] t time at which to evaluate residual \f$t^n\f$
 * \param[out] ss_flux steady-state rhs vector \f$\mathrm{q}^n\f$
 */
template <int dim>
void ShallowWater<dim>::compute_ss_rhs(const double &, Vector<double> & ss_rhs)
{
  ss_rhs = 0;
}

/**
 * \brief Computes the inviscid fluxes required to be evaluated in cell and face
 *        integrations.
 *
 * \param[in] height  vector of height values
 * \param[in] momentum vector of momentum values
 * \param[out] height_flux  vector of height inviscid flux values
 * \param[out] momentum_flux vector of momentum inviscid flux values
*/
template <int dim>
void ShallowWater<dim>::compute_inviscid_fluxes(
  const std::vector<double> & height,
  const std::vector<Tensor<1, dim>> & momentum,
  std::vector<Tensor<1, dim>> & height_flux,
  std::vector<Tensor<2, dim>> & momentum_flux) const
{
  // get number of vector elements
  const unsigned int n = height.size();

  // identity tensor
  SymmetricTensor<2, dim> identity_tensor_sym = unit_symmetric_tensor<dim>();
  Tensor<2, dim> identity_tensor(identity_tensor_sym);

  // compute auxiliary quantities
  std::vector<Tensor<1, dim>> velocity = compute_velocity(height, momentum);

  // loop over vector elements
  for (unsigned int q = 0; q < n; ++q)
  {
    // compute height inviscid flux
    height_flux[q] = momentum[q];

    // compute momentum inviscid flux
    Tensor<2, dim> velocity_times_momentum =
      outer_product(velocity[q], momentum[q]);
    momentum_flux[q] = velocity_times_momentum +
      0.5 * gravity * std::pow(height[q], 2) * identity_tensor;
  }
}

/**
 * \brief Computes a vector of velocity values.
 *
 * \param[in] height   vector of height values
 * \param[in] momentum vector of momentum values
 *
 * \return vector of velocity values
 */
template <int dim>
std::vector<Tensor<1, dim>> ShallowWater<dim>::compute_velocity(
  const std::vector<double> & height,
  const std::vector<Tensor<1, dim>> & momentum) const
{
  const unsigned int n = height.size();

  std::vector<Tensor<1, dim>> velocity(n);
  for (unsigned int q = 0; q < n; ++q)
    velocity[q] = momentum[q] / height[q];

  return velocity;
}

/**
 * \brief Computes a vector of speed values.
 *
 * \param[in] height   vector of height values
 * \param[in] momentum vector of momentum values
 *
 * \return vector of speed values
 */
template <int dim>
std::vector<double> ShallowWater<dim>::compute_speed(
  const std::vector<double> & height,
  const std::vector<Tensor<1, dim>> & momentum) const
{
  const unsigned int n = height.size();

  std::vector<double> speed(n);
  for (unsigned int q = 0; q < n; ++q)
    speed[q] = momentum[q].norm() / height[q];

  return speed;
}

/**
 * \brief Computes a vector of sound speed values.
 *
 * \param[in] height vector of height values
 *
 * \return vector of sound speed values
 */
template <int dim>
std::vector<double> ShallowWater<dim>::compute_sound_speed(
  const std::vector<double> & height) const
{
  const unsigned int n = height.size();

  std::vector<double> sound_speed(n);
  for (unsigned int q = 0; q < n; ++q)
    sound_speed[q] = std::sqrt(gravity * height[q]);

  return sound_speed;
}

/**
 * \brief Computes the maximum flux speed \f$\lambda\f$ at each quadrature point
 *        in domain and finds the max in each cell and the max in the
 *        entire domain.
 *
 * The eigenvalues of the shallow water equations are:
 * - \f$\lambda_1 = \|\mathbf{v}\|\f$,
 * - \f$\lambda_2 = \|\mathbf{v}\| + a\f$,
 * - \f$\lambda_3 = \|\mathbf{v}\| - a\f$,
 *
 * where \f$a\f$ is the sound speed:
 * \f[
 *   a = \sqrt{gh} \,.
 * \f]
 * The maximum wave speed is thus the second eigenvalue:
 * \f[
 *   \lambda_{max} = \|\mathbf{v}\| + a \,.
 * \f]
 */
template <int dim>
void ShallowWater<dim>::update_flux_speeds()
{
  FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);

  // reset max flux speed
  this->max_flux_speed = 0.0;

  // loop over cells
  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    fe_values.reinit(cell);

    // get solution values for height and momentum
    std::vector<double> height(this->n_q_points_cell);
    std::vector<Tensor<1, dim>> momentum(this->n_q_points_cell);
    fe_values[height_extractor].get_function_values(this->new_solution, height);
    fe_values[momentum_extractor].get_function_values(this->new_solution,
                                                      momentum);

    // compute velocity values
    std::vector<Tensor<1, dim>> velocity = compute_velocity(height, momentum);

    // determine max flux speed in cell
    this->max_flux_speed_cell[cell] = 0.0;
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
      // compute flux speed
      double flux_speed = velocity[q].norm() + std::sqrt(gravity * height[q]);

      // compare with current maximum flux speed in cell
      this->max_flux_speed_cell[cell] =
        std::max(this->max_flux_speed_cell[cell], flux_speed);
    }

    // determine max flux speed in domain
    this->max_flux_speed =
      std::max(this->max_flux_speed, this->max_flux_speed_cell[cell]);
  }
}

/**
 * \brief Returns a pointer to a shallow water auxiliary post-processor object
 *
 * \return pointer to a shallow water auxiliary post-processor object
 */
template <int dim>
std::shared_ptr<DataPostprocessor<dim>> ShallowWater<
  dim>::create_auxiliary_postprocessor() const
{
  return std::make_shared<ShallowWaterPostProcessor<dim>>(bathymetry_function);
}

/**
 * \brief Creates a viscosity multiplier object that uses the local Froude
 *        number as the multiplier.
 *
 * \return pointer to created viscosity multiplier
 */
template <int dim>
std::shared_ptr<ViscosityMultiplier<dim>> ShallowWater<
  dim>::create_viscosity_multiplier() const
{
  std::shared_ptr<ViscosityMultiplier<dim>> viscosity_multiplier;
  if (sw_parameters.multiply_low_order_viscosity_by_froude)
  {
    viscosity_multiplier =
      std::make_shared<ShallowWaterViscosityMultiplier<dim>>(gravity);
  }
  else
  {
    viscosity_multiplier = std::make_shared<ViscosityMultiplier<dim>>();
  }

  return viscosity_multiplier;
}
