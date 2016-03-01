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
ShallowWater<dim>::ShallowWater(const ShallowWaterRunParameters<dim> & params)
  : ConservationLaw<dim>(params, dim + 1, true),
    sw_parameters(params),
    problem_parameters(params.problem_name, false),
    height_extractor(0),
    momentum_extractor(1),
    fe_bathymetry(params.degree),
    dof_handler_bathymetry(this->triangulation)
{
  // shallow water equations cannot be 3-D
  Assert(dim < 3, ExcImpossibleInDim(dim));

  // point base class problem parameters to derived class problem parameters
  this->problem_base_parameters = &problem_parameters;
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
  // get path of source directory from #define created by CMake
  std::stringstream source_path_ss;
  source_path_ss << SOURCE_PATH;
  std::string source_path;
  source_path_ss >> source_path;

  // create problem parameters file name and determine if it exists
  std::string problem_parameters_file =
    source_path + "/problems/shallowwater/" + sw_parameters.problem_name;

  // get and process the problem parameters
  problem_parameters.get_and_process_parameters(problem_parameters_file,
                                                this->triangulation,
                                                this->fe,
                                                this->face_quadrature);
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
                                               problem_parameters.gravity,
                                               bathymetry_vector,
                                               problem_parameters.domain_volume,
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
    std::make_shared<ShallowWaterStarState<dim>>(problem_parameters.gravity,
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
  auto max_wave_speed = std::make_shared<ShallowWaterMaxWaveSpeed<dim>>(
    this->star_state, problem_parameters.gravity);

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
  VectorTools::interpolate(dof_handler_bathymetry,
                           *(problem_parameters.bathymetry_function),
                           bathymetry_vector);
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
            fe_values[momentum_extractor].value(i, q) *
              problem_parameters.gravity * height[q] * bathymetry_gradient[q]) *
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
      0.5 * problem_parameters.gravity * std::pow(height[q], 2) * identity_tensor;
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
    sound_speed[q] = std::sqrt(problem_parameters.gravity * height[q]);

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
      double flux_speed =
        velocity[q].norm() + std::sqrt(problem_parameters.gravity * height[q]);

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
  return std::make_shared<ShallowWaterPostProcessor<dim>>(
    problem_parameters.bathymetry_function);
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
    viscosity_multiplier = std::make_shared<ShallowWaterViscosityMultiplier<dim>>(
      problem_parameters.gravity);
  }
  else
  {
    viscosity_multiplier = std::make_shared<ViscosityMultiplier<dim>>();
  }

  return viscosity_multiplier;
}

/**
 * \brief Creates an FCT object and returns the pointer.
 *
 * \return FCT object pointer
 */
template <int dim>
std::shared_ptr<FCT<dim>> ShallowWater<dim>::create_fct() const
{
  // determine if star states are to be used in FCT bounds
  const bool use_star_states_in_fct_bounds =
    this->are_star_states && this->parameters.use_star_states_in_fct_bounds;

  auto fct =
    std::make_shared<ShallowWaterFCT<dim>>(this->parameters,
                                           this->dof_handler,
                                           this->triangulation,
                                           this->lumped_mass_matrix,
                                           this->consistent_mass_matrix,
                                           this->star_state,
                                           *(this->linear_solver),
                                           this->unconstrained_sparsity_pattern,
                                           this->dirichlet_dof_indices,
                                           this->n_components,
                                           this->dofs_per_cell,
                                           this->component_names,
                                           use_star_states_in_fct_bounds,
                                           problem_parameters.gravity);

  return fct;
}
