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
  : ConservationLaw<dim>(params),
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
void ShallowWater<dim>::define_problem()
{
  // determine problem name
  this->problem_name = sw_parameters.problem_name;

  // create problem parameters file name and determine if it exists
  std::string problem_parameters_file =
    "problems/shallow_water/" + this->problem_name;
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
  this->n_boundaries = 1;
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
    std::shared_ptr<DirichletBoundaryConditions<dim>>
      derived_boundary_conditions =
        std::make_shared<DirichletBoundaryConditions<dim>>(this->fe,
                                                           this->face_quadrature);
    this->boundary_conditions = derived_boundary_conditions;

    // option to use exact solution as Dirichlet BC
    this->use_exact_solution_as_dirichlet_bc = false;

    // get Dirichlet function strings
    if (!this->use_exact_solution_as_dirichlet_bc)
    {
      this->dirichlet_function_strings.resize(this->n_boundaries);
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
 * \brief Interpolates the bathymetry FE vector from its function.
 */
template <int dim>
void ShallowWater<dim>::perform_additional_setup()
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

    // add to global mass matrix with contraints
    this->constraints.distribute_local_to_global(
      local_mass, local_dof_indices, this->lumped_mass_matrix);
  }
}

/**
 * \brief Computes the steady-state residual.
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
 *  \param[in] dt time step size \f$\Delta t\f$
 *  \param[out] r steady-state residual \f$\mathbf{r}\f$
 */
template <int dim>
void ShallowWater<dim>::compute_ss_residual(const double & dt, Vector<double> & f)
{
  // reset vector
  f = 0.0;

  FEValues<dim> fe_values(this->fe,
                          this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  FEValues<dim> fe_values_bathymetry(
    fe_bathymetry, this->cell_quadrature, update_gradients);

  Vector<double> cell_residual(this->dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);

  // loop over cells
unsigned int i_cell = 0;
  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end(),
       cell_bathymetry = dof_handler_bathymetry.begin_active();
  for (; cell != endc; ++cell, ++cell_bathymetry, ++i_cell)
  {
std::cout << "visc[" << i_cell << "] = " << this->viscosity[cell] << std::endl;
    // reset cell residual
    cell_residual = 0;

    // reinitialize fe values for cell
    fe_values.reinit(cell);
    fe_values_bathymetry.reinit(cell_bathymetry);

    // get current solution values
    std::vector<double> height(this->n_q_points_cell);
    std::vector<Tensor<1, dim>> momentum(this->n_q_points_cell);
    fe_values[height_extractor].get_function_values(this->new_solution, height);
    fe_values[momentum_extractor].get_function_values(this->new_solution,
                                                      momentum);

    // get solution gradients
    std::vector<Tensor<1, dim>> height_gradient(this->n_q_points_cell);
    std::vector<Tensor<2, dim>> momentum_gradient(this->n_q_points_cell);
    fe_values[height_extractor].get_function_gradients(this->new_solution,
                                                       height_gradient);
    fe_values[momentum_extractor].get_function_gradients(this->new_solution,
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

    // compute viscous fluxes
    std::vector<Tensor<1, dim>> height_viscous_flux(this->n_q_points_cell);
    std::vector<Tensor<2, dim>> momentum_viscous_flux(this->n_q_points_cell);
    compute_viscous_fluxes(this->viscosity[cell],
                           height_gradient,
                           momentum_gradient,
                           bathymetry_gradient,
                           height_viscous_flux,
                           momentum_viscous_flux);

    // get quadrature points on cell
    std::vector<Point<dim>> points(this->n_q_points_cell);
    points = fe_values.get_quadrature_points();

    // loop over quadrature points
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        cell_residual(i) +=
          (
            // height inviscid flux
            fe_values[height_extractor].gradient(i, q) *
              (height_inviscid_flux[q] + height_viscous_flux[q])
            // momentum inviscid flux
            +
            double_contract<0, 0, 1, 1>(
              fe_values[momentum_extractor].gradient(i, q),
              momentum_inviscid_flux[q] + momentum_viscous_flux[q])
            // bathymetry source term
            -
            fe_values[momentum_extractor].value(i, q) * gravity * height[q] *
              bathymetry_gradient[q]) *
          fe_values.JxW(q);
      }
    }

    // apply boundary conditions
    this->boundary_conditions->apply(
      cell, fe_values, this->new_solution, dt, cell_residual);

    // aggregate local residual into global residual
    cell->get_dof_indices(local_dof_indices);
    this->constraints.distribute_local_to_global(
      cell_residual, local_dof_indices, f);
  } // end cell loop
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
 * \brief Computes the viscous fluxes.
 */
template <int dim>
void ShallowWater<dim>::compute_viscous_fluxes(
  const double & viscosity,
  const std::vector<Tensor<1, dim>> & height_gradient,
  const std::vector<Tensor<2, dim>> & momentum_gradient,
  const std::vector<Tensor<1, dim>> & bathymetry_gradient,
  std::vector<Tensor<1, dim>> & height_viscous_flux,
  std::vector<Tensor<2, dim>> & momentum_viscous_flux) const
{
  // get number of vector elements
  const unsigned int n = height_gradient.size();

  // loop over vector elements
  for (unsigned int q = 0; q < n; ++q)
  {
    // density viscous flux
    height_viscous_flux[q] =
      -viscosity * (height_gradient[q] + 0 * bathymetry_gradient[q]);

    // momentum viscous flux
    momentum_viscous_flux[q] = -viscosity * momentum_gradient[q];
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

  // loop over cells to compute first order viscosity at each quadrature point
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
 * \brief Computes entropy \f$\eta\f$ at each quadrature point on cell or face.
 *
 * For the shallow water equations, the entropy is defined as
 * \f[
 *   \eta(\mathbf{u}) = \frac{1}{2}\frac{\mathbf{q}\cdot\mathbf{q}}{h}
 *   + \frac{1}{2}g h^2
 * \f]
 */
template <int dim>
void ShallowWater<dim>::compute_entropy(const Vector<double> & solution,
                                        const FEValuesBase<dim> & fe_values,
                                        Vector<double> & entropy) const
{
  // get number of quadrature points
  const unsigned int n = entropy.size();

  // get height and momentum
  std::vector<double> height(n);
  std::vector<Tensor<1, dim>> momentum(n);
  fe_values[height_extractor].get_function_values(solution, height);
  fe_values[momentum_extractor].get_function_values(solution, momentum);

  // compute entropy
  for (unsigned int q = 0; q < n; ++q)
    entropy(q) = 0.5 * momentum[q] * momentum[q] / height[q] +
      0.5 * gravity * height[q] * height[q];
}

/**
 * \brief Computes low-order viscosity for each cell.
 *
 * The low-order viscosity is computed as
 * \f[
 *   \nu_K^L = c_{max} h_K \lambda_{K,max} \,,
 * \f]
 * where \f$\lambda_{K,max}\f$ is the maximum flux speed on cell \f$K\f$,
 * or if the user specifies, the low-order viscosity definition above is
 * multiplied by the local Froude number:
 * \f[
 *   \nu_K^L = c_{max} \mbox{Fr}_K h_K \lambda_{K,max} \,,
 * \f]
 * where the Froude number is taken to be the max over all quadrature points:
 * \f[
 *   \mbox{Fr}_K = \max\limits_q \frac{u_q}{a_q} \,.
 * \f]
 *
 * \param[in] using_low_order_scheme flag that the low-order viscosities are
 *            are used in a low-order scheme as opposed to an entropy
 *            viscosity scheme
 */
template <int dim>
void ShallowWater<dim>::update_old_low_order_viscosity(
  const bool & using_low_order_scheme)
{
  const double c_max = this->parameters.first_order_viscosity_coef;

  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();

  if (using_low_order_scheme &&
      sw_parameters.multiply_low_order_viscosity_by_froude)
  {
    FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);

    // compute low-order viscosity for each cell
    for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      // extract height and momentum
      std::vector<double> height(this->n_q_points_cell);
      fe_values[height_extractor].get_function_values(this->new_solution, height);
      std::vector<Tensor<1, dim>> momentum(this->n_q_points_cell);
      fe_values[momentum_extractor].get_function_values(this->new_solution,
                                                        momentum);

      // compute fluid speed and sound speed
      std::vector<double> speed = compute_speed(height, momentum);
      std::vector<double> sound_speed = compute_sound_speed(height);

      // compute Froude number for cell by taking the maximum over
      // quadrature points
      double froude = 0.0;
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
        froude = std::max(froude, speed[q] / sound_speed[q]);

      this->first_order_viscosity[cell] =
        std::abs(c_max * froude * this->cell_diameter[cell] *
                 this->max_flux_speed_cell[cell]);
    }
  }
  else
  {
    // compute low-order viscosity for each cell
    for (; cell != endc; ++cell)
    {
      this->first_order_viscosity[cell] = std::abs(
        c_max * this->cell_diameter[cell] * this->max_flux_speed_cell[cell]);
    }
  }
}

/**
 * \brief Computes entropy viscosity for each cell.
 *
 * \param[in] dt time step size
 */
template <int dim>
void ShallowWater<dim>::update_entropy_viscosities(const double & dt)
{
  // compute normalization constant for entropy viscosity
  double entropy_average = 1.0e15;
  if (sw_parameters.entropy_normalization == "average")
    entropy_average = this->compute_average_entropy(this->new_solution);

  // get tuning parameters
  const double entropy_residual_coefficient =
    this->parameters.entropy_residual_coef;
  const double jump_coefficient = this->parameters.entropy_jump_coef;

  // FE values for entropy flux
  ShallowWaterEntropyFluxFEValuesCell<dim> entropy_flux_fe_values_cell(
    this->dof_handler,
    this->triangulation,
    this->cell_quadrature,
    this->new_solution,
    bathymetry_vector,
    gravity);
  ShallowWaterEntropyFluxFEValuesFace<dim> entropy_flux_fe_values_face(
    this->dof_handler,
    this->triangulation,
    this->face_quadrature,
    this->new_solution,
    bathymetry_vector,
    gravity);

  // compute entropy viscosity for each cell
unsigned int i_cell = 0;
  Cell cell = this->dof_handler.begin_active(), endc = this->dof_handler.end();
  for (; cell != endc; ++cell, ++i_cell)
  {
    // reinitialize entropy flux FE values for cell (face values will need be
    // reinitialized in compute_max_entropy_jump())
    entropy_flux_fe_values_cell.reinit(cell);

    // compute entropy normalization
    std::vector<double> entropy_normalization(this->n_q_points_cell);
    if (sw_parameters.entropy_normalization == "average")
      entropy_normalization = this->compute_entropy_normalization(
        this->new_solution, entropy_average, cell);
    else if (sw_parameters.entropy_normalization == "local")
      entropy_normalization =
        compute_local_entropy_normalization(this->new_solution, cell);
    else if (sw_parameters.entropy_normalization == "constant")
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
        entropy_normalization[q] =
          sw_parameters.constant_entropy_normalization_coefficient;
    else
      ExcNotImplemented();

    // compute entropy residual at each quadrature point on cell
    const std::vector<double> entropy_residual =
      this->compute_entropy_residual(this->new_solution,
                                     this->old_solution,
                                     entropy_flux_fe_values_cell,
                                     dt,
                                     cell);

    // compute max entropy flux jump
    const double max_entropy_jump =
      compute_max_entropy_jump(entropy_flux_fe_values_face, cell);
std::cout << "jump[" << i_cell << "] = " << max_entropy_jump << std::endl;

    // compute entropy viscosity
    double h2 = std::pow(this->cell_diameter[cell], 2);
    this->entropy_viscosity[cell] = 0.0;
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      this->entropy_viscosity[cell] = std::max(
        this->entropy_viscosity[cell],
        h2 * (entropy_residual_coefficient * std::abs(entropy_residual[q]) +
              jump_coefficient * max_entropy_jump) /
          entropy_normalization[q]);
  }
}

/**
 * \brief Computes the local entropy viscosity normalization coefficient
 *        \f$gh^2\f$ for each quadrature point in a cell.
 *
 * The local entropy viscosity normalization coefficient is computed as
 *
 * \f[
 *   c^{\mbox{normalization}}_q = g h_q^2
 * \f]
 *
 * \param[in] solution solution vector
 * \param[in] cell cell iterator
 *
 * \return vector of entropy viscosity normalization coefficient for each
 *         quadrature point in cell
 */
template <int dim>
std::vector<double> ShallowWater<dim>::compute_local_entropy_normalization(
  const Vector<double> & solution, const Cell & cell) const
{
  // get height values
  FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);
  fe_values.reinit(cell);
  std::vector<double> height(this->n_q_points_cell);
  fe_values[height_extractor].get_function_values(solution, height);

  // compute normalization at each quadrature point
  std::vector<double> normalization(this->n_q_points_cell);
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    normalization[q] = gravity * std::pow(height[q], 2);

  return normalization;
}

/**
 * \brief Computes the entropy residual at each quadrature point in a cell.
 *
 * \param[in] new_solution new solution
 * \param[in] old_solution old solution
 * \param[in] entropy_flux_fe_values FE values for entropy flux
 * \param[in] dt time step size
 * \param[in] cell cell iterator
 *
 * \return vector of entropy residual for each quadrature point in cell
 */
template <int dim>
std::vector<double> ShallowWater<dim>::compute_entropy_residual(
  const Vector<double> & new_solution,
  const Vector<double> & old_solution,
  const ShallowWaterEntropyFluxFEValuesCell<dim> & entropy_flux_fe_values,
  const double & dt,
  const Cell & cell) const
{
  // FE values
  FEValues<dim> fe_values(
    this->fe, this->cell_quadrature, update_values | update_gradients);
  fe_values.reinit(cell);

  Vector<double> entropy_new(this->n_q_points_cell);
  Vector<double> entropy_old(this->n_q_points_cell);
  std::vector<double> entropy_residual(this->n_q_points_cell);

  // compute entropy of current and old solutions
  compute_entropy(new_solution, fe_values, entropy_new);
  compute_entropy(old_solution, fe_values, entropy_old);
  std::vector<double> divergence_entropy_flux =
    entropy_flux_fe_values.get_function_divergences();

  // compute entropy residual at each quadrature point on cell
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
  {
    // compute entropy residual
    double dsdt = (entropy_new[q] - entropy_old[q]) / dt;
    entropy_residual[q] = dsdt + divergence_entropy_flux[q];
  }

  return entropy_residual;
}

/**
 * \brief Computes the max entropy jump in a cell.
 *
 * \param[in] fe_values entropy flux FE values
 * \param[in] cell cell iterator
 *
 * \return max entropy jump in cell
 */
template <int dim>
double ShallowWater<dim>::compute_max_entropy_jump(
  ShallowWaterEntropyFluxFEValuesFace<dim> & fe_values, const Cell & cell) const
{
  std::vector<Tensor<2, dim>> entropy_flux_gradients_this_cell(
    this->n_q_points_face);
  std::vector<Tensor<2, dim>> entropy_flux_gradients_neighbor_cell(
    this->n_q_points_face);
  std::vector<Tensor<1, dim>> normal_vectors(this->n_q_points_face);

  double max_jump_in_cell = 0.0;

  // loop over faces
  for (unsigned int iface = 0; iface < this->faces_per_cell; ++iface)
  {
    // determine if face is interior
    typename DoFHandler<dim>::face_iterator face = cell->face(iface);
    if (face->at_boundary() == false)
    {
      // reinitialize FE values
      fe_values.reinit(cell, iface);

      // get gradients
      fe_values.get_function_gradients(entropy_flux_gradients_this_cell);

      // get normal vectors
      normal_vectors = fe_values.get_normal_vectors();

      // get iterator to neighboring cell and face index of this face
      Assert(cell->neighbor(iface).state() == IteratorState::valid,
             ExcInternalError());
      Cell neighbor = cell->neighbor(iface);
      const unsigned int ineighbor = cell->neighbor_of_neighbor(iface);
      Assert(ineighbor < this->faces_per_cell, ExcInternalError());

      // get gradients from neighboring cell
      fe_values.reinit(neighbor, ineighbor);
      fe_values.get_function_gradients(entropy_flux_gradients_neighbor_cell);

      // loop over face quadrature points to determine max jump on face
      double max_jump_on_face = 0.0;
      for (unsigned int q = 0; q < this->n_q_points_face; ++q)
      {
        // compute difference in gradients across face
        Tensor<2, dim> entropy_flux_gradient_jump =
          entropy_flux_gradients_this_cell[q] -
          entropy_flux_gradients_neighbor_cell[q];
        double jump_on_face = std::abs(entropy_flux_gradient_jump *
                                       normal_vectors[q] * normal_vectors[q]);
        max_jump_on_face = std::max(max_jump_on_face, jump_on_face);
      }

      // update max jump in cell
      max_jump_in_cell = std::max(max_jump_in_cell, max_jump_on_face);
    }
  }

  return max_jump_in_cell;
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
