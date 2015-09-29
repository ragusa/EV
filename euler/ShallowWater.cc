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
  switch (sw_parameters.problem_id)
  {
    case 0: // 1-D dam break over flat bottom
    {
      Assert(dim == 1, ExcImpossibleInDim(dim));

      // name of problem
      this->problem_name = "dam_break_flat";

      // domain
      const double domain_start = -5.0;
      const double domain_width = 10.0;
      this->domain_volume = std::pow(domain_width, dim);
      GridGenerator::hyper_cube(
        this->triangulation, domain_start, domain_start + domain_width);

      // only 1 type of BC: zero Dirichlet; leave boundary indicators as zero
      this->n_boundaries = 1;
      typename Triangulation<dim>::cell_iterator cell =
        this->triangulation.begin();
      typename Triangulation<dim>::cell_iterator endc = this->triangulation.end();
      for (; cell != endc; ++cell)
        for (unsigned int face = 0; face < this->faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_id(0);
      this->boundary_types.resize(this->n_boundaries);
      this->boundary_types[0].resize(this->n_components);
      this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet;
      this->boundary_types[0][1] = ConservationLaw<dim>::dirichlet;
      this->dirichlet_function_strings.resize(this->n_boundaries);
      this->dirichlet_function_strings[0].resize(this->n_components);
      this->dirichlet_function_strings[0][0] = "if(x<0,3,1)";
      this->dirichlet_function_strings[0][1] = "0";
      this->use_exact_solution_as_BC = false;

      // initial conditions
      this->initial_conditions_strings[0] = "if(x<0,3,1)";
      this->initial_conditions_strings[1] = "0";

      // problem parameters
      const double h_left = 3.0;
      const double u_left = 0.0;
      const double h_right = 1.0;
      const double u_right = 0.0;
      gravity = 1.0;
      const double x_interface = 0.0;

      // exact solution
      this->has_exact_solution = true;
      // create and initialize Riemann solver for exact solution
      std::shared_ptr<ShallowWaterRiemannSolver<dim>>
        exact_solution_function_derived =
          std::make_shared<ShallowWaterRiemannSolver<dim>>(
            h_left, u_left, h_right, u_right, gravity, x_interface);
      // point base class pointer to derived class function object
      this->exact_solution_function = exact_solution_function_derived;

      // initialize bathymetry function
      std::shared_ptr<FunctionParser<dim>> bathymetry_function_derived =
        std::make_shared<FunctionParser<dim>>();
      std::map<std::string, double> constants;
      std::string bathymetry_string = "0";
      bathymetry_function_derived->initialize(
        "x", bathymetry_string, this->constants, false);
      bathymetry_function = bathymetry_function_derived;

      // initialize bathymetry gradient function
/*
      std::shared_ptr<FunctionParser<dim>> bathymetry_gradient_function_derived =
        std::make_shared<FunctionParser<dim>>(dim);
      std::string bathymetry_gradient_string = "0";
      bathymetry_gradient_function_derived->initialize(
        "x", bathymetry_gradient_string, this->constants, false);
      bathymetry_gradient_function = bathymetry_gradient_function_derived;
*/

      // default end time
      this->has_default_end_time = true;
      this->default_end_time = 2.0;

      break;
    }
    case 1: // 1-D dam break over rectangular bump
    {
      Assert(dim == 1, ExcImpossibleInDim(dim));

      // NOTE: need to ensure that adaptive mesh refinement is not used for
      // this problem because bump edges must be coincident with mesh edges;
      // otherwise, the gradient of the discontinuity would be evaluated

      // name of problem
      this->problem_name = "dam_break_bump";

      // domain
      const double domain_start = 0.0;
      const double domain_width = 1500.0;
      this->domain_volume = std::pow(domain_width, dim);
      GridGenerator::hyper_cube(
        this->triangulation, domain_start, domain_start + domain_width);

      // constants
      const double midpoint = domain_start + 0.5 * domain_width;
      const double bump_width = 375.0;
      const double bump_left_nominal = midpoint - 0.5 * bump_width;
      const double bump_right_nominal = midpoint + 0.5 * bump_width;

      // determine actual bump left and right positions; they must be coincident
      // with mesh edges
      const unsigned int n_cells =
        std::pow(2, this->parameters.initial_refinement_level);
      const double dx = domain_width / n_cells;
      const double bump_left =
        std::round((bump_left_nominal - domain_start) / dx) * dx;
      const double bump_right =
        std::round((bump_right_nominal - domain_start) / dx) * dx;

      // constants for function parser
      this->constants["midpoint"] = midpoint;     // midpoint
      this->constants["bump_width"] = bump_width; // width of bump
      this->constants["bump_height"] = 8.0;       // height of bump
      this->constants["bump_left"] = bump_left;   // left edge of bump
      this->constants["bump_right"] = bump_right; // right edge of bump
      this->constants["h_left"] = 20.0;
      this->constants["h_right"] = 15.0;

      // only 1 type of BC: zero Dirichlet; leave boundary indicators as zero
      this->n_boundaries = 1;
      typename Triangulation<dim>::cell_iterator cell =
        this->triangulation.begin();
      typename Triangulation<dim>::cell_iterator endc = this->triangulation.end();
      for (; cell != endc; ++cell)
        for (unsigned int face = 0; face < this->faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_id(0);
      this->boundary_types.resize(this->n_boundaries);
      this->boundary_types[0].resize(this->n_components);
      this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet;
      this->boundary_types[0][1] = ConservationLaw<dim>::dirichlet;
      this->dirichlet_function_strings.resize(this->n_boundaries);
      this->dirichlet_function_strings[0].resize(this->n_components);
      this->dirichlet_function_strings[0][0] = "if(x<midpoint,h_left,h_right)";
      this->dirichlet_function_strings[0][1] = "0";
      this->use_exact_solution_as_BC = false;

      // initial conditions
      this->initial_conditions_strings[0] =
        "if(x<bump_left,h_left,"
        "if(x<midpoint,h_left-bump_height,"
        "if(x<=bump_right,h_right-bump_height,h_right)))";
      this->initial_conditions_strings[1] = "0";

      // problem parameters
      gravity = 9.812;

      // exact solution
      this->has_exact_solution = false;

      // initialize bathymetry function
      std::shared_ptr<FunctionParser<dim>> bathymetry_function_derived =
        std::make_shared<FunctionParser<dim>>();
      std::string bathymetry_string = "if(x<bump_left,0,"
                                      "if(x<=bump_right,bump_height,0))";
      bathymetry_function_derived->initialize(
        "x", bathymetry_string, this->constants, false);
      bathymetry_function = bathymetry_function_derived;

      // initialize bathymetry gradient function
/*
      std::shared_ptr<FunctionParser<dim>> bathymetry_gradient_function_derived =
        std::make_shared<FunctionParser<dim>>(dim);
      std::string bathymetry_gradient_string = "0";
      bathymetry_gradient_function_derived->initialize(
        "x", bathymetry_gradient_string, this->constants, false);
      bathymetry_gradient_function = bathymetry_gradient_function_derived;
*/

      // default end time
      this->has_default_end_time = true;
      this->default_end_time = 15.0;

      break;
    }
    case 2: // 1-D lake at rest
    {
      Assert(dim == 1, ExcImpossibleInDim(dim));

      // name of problem
      this->problem_name = "lake_at_rest";

      // domain
      const double domain_start = 0.0;
      const double domain_width = 20.0;
      this->domain_volume = std::pow(domain_width, dim);
      GridGenerator::hyper_cube(
        this->triangulation, domain_start, domain_start + domain_width);

      // only 1 type of BC: zero Dirichlet; leave boundary indicators as zero
      this->n_boundaries = 1;
      typename Triangulation<dim>::cell_iterator cell =
        this->triangulation.begin();
      typename Triangulation<dim>::cell_iterator endc = this->triangulation.end();
      for (; cell != endc; ++cell)
        for (unsigned int face = 0; face < this->faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_id(0);
      this->boundary_types.resize(this->n_boundaries);
      this->boundary_types[0].resize(this->n_components);
      this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet;
      this->boundary_types[0][1] = ConservationLaw<dim>::dirichlet;
      this->dirichlet_function_strings.resize(this->n_boundaries);
      this->dirichlet_function_strings[0].resize(this->n_components);
      this->dirichlet_function_strings[0][0] = "1"; // height
      this->dirichlet_function_strings[0][1] = "0"; // momentum
      this->use_exact_solution_as_BC = false;

      // initial conditions
      this->initial_conditions_strings[0] = "if(abs(x-10)<2,1-(4-(x-10)^2)/20,1)";
      this->initial_conditions_strings[1] = "0";

      // problem parameters
      gravity = 9.812;

      // exact solution
      this->has_exact_solution = true;
      this->exact_solution_strings[0] = "if(abs(x-10)<2,1-(4-(x-10)^2)/20,1)";
      this->exact_solution_strings[1] = "0";

      // create and initialize function parser for exact solution
      std::shared_ptr<FunctionParser<dim>> exact_solution_function_derived =
        std::make_shared<FunctionParser<dim>>(this->parameters.n_components);
      exact_solution_function_derived->initialize(
        "x", this->exact_solution_strings, this->constants, false);
      this->exact_solution_function = exact_solution_function_derived;

      // initialize bathymetry function
      std::shared_ptr<FunctionParser<dim>> bathymetry_function_derived =
        std::make_shared<FunctionParser<dim>>();
      std::string bathymetry_string = "if(abs(x-10)<2,(4-(x-10)^2)/20,0)";
      bathymetry_function_derived->initialize(
        "x", bathymetry_string, this->constants, false);
      bathymetry_function = bathymetry_function_derived;

      // initialize bathymetry gradient function
/*
      std::shared_ptr<FunctionParser<dim>> bathymetry_gradient_function_derived =
        std::make_shared<FunctionParser<dim>>(dim);
      std::string bathymetry_gradient_string = "if(abs(x-10)<2,1-x/10,0)";
      bathymetry_gradient_function_derived->initialize(
        "x", bathymetry_gradient_string, this->constants, false);
      bathymetry_gradient_function = bathymetry_gradient_function_derived;
*/

      // default end time
      this->has_default_end_time = true;
      this->default_end_time = 10.0;

      break;
    }
    case 3: // 1-D lake at rest perturbed
    {
      Assert(dim == 1, ExcImpossibleInDim(dim));

      // name of problem
      this->problem_name = "lake_at_rest_perturbed";

      // domain
      const double domain_start = 0.0;
      const double domain_width = 20.0;
      this->domain_volume = std::pow(domain_width, dim);
      GridGenerator::hyper_cube(
        this->triangulation, domain_start, domain_start + domain_width);

      // only 1 type of BC: zero Dirichlet; leave boundary indicators as zero
      this->n_boundaries = 1;
      typename Triangulation<dim>::cell_iterator cell =
        this->triangulation.begin();
      typename Triangulation<dim>::cell_iterator endc = this->triangulation.end();
      for (; cell != endc; ++cell)
        for (unsigned int face = 0; face < this->faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_id(0);
      this->boundary_types.resize(this->n_boundaries);
      this->boundary_types[0].resize(this->n_components);
      this->boundary_types[0][0] = ConservationLaw<dim>::dirichlet;
      this->boundary_types[0][1] = ConservationLaw<dim>::dirichlet;
      this->dirichlet_function_strings.resize(this->n_boundaries);
      this->dirichlet_function_strings[0].resize(this->n_components);
      this->dirichlet_function_strings[0][0] = "1"; // height
      this->dirichlet_function_strings[0][1] = "0"; // momentum
      this->use_exact_solution_as_BC = false;

      // initial conditions
      this->initial_conditions_strings[0] = "if(abs(x-6)<0.25,1.01,"
        "if(abs(x-10)<2,1-(4-(x-10)^2)/20,1))";
      this->initial_conditions_strings[1] = "0";

      // problem parameters
      gravity = 9.812;

      // exact solution
      this->has_exact_solution = true;
      this->exact_solution_strings[0] = "if(abs(x-10)<2,1-(4-(x-10)^2)/20,1)";
      this->exact_solution_strings[1] = "0";

      // create and initialize function parser for exact solution
      std::shared_ptr<FunctionParser<dim>> exact_solution_function_derived =
        std::make_shared<FunctionParser<dim>>(this->parameters.n_components);
      exact_solution_function_derived->initialize(
        "x", this->exact_solution_strings, this->constants, false);
      this->exact_solution_function = exact_solution_function_derived;

      // initialize bathymetry function
      std::shared_ptr<FunctionParser<dim>> bathymetry_function_derived =
        std::make_shared<FunctionParser<dim>>();
      std::string bathymetry_string = "if(abs(x-10)<2,(4-(x-10)^2)/20,0)";
      bathymetry_function_derived->initialize(
        "x", bathymetry_string, this->constants, false);
      bathymetry_function = bathymetry_function_derived;

      // default end time
      this->has_default_end_time = true;
      this->default_end_time = 1.5;

      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }
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

  cell_iterator cell = this->dof_handler.begin_active(),
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
 *  \param[out] r steady-state residual \f$\mathbf{r}\f$
 */
template <int dim>
void ShallowWater<dim>::compute_ss_residual(Vector<double> & f)
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
  cell_iterator cell = this->dof_handler.begin_active(),
                endc = this->dof_handler.end(),
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
                           height_viscous_flux,
                           momentum_viscous_flux);

    // compute gradients of bathymetry function
    std::vector<Tensor<1, dim>> bathymetry_gradient(this->n_q_points_cell);
    fe_values_bathymetry.get_function_gradients(bathymetry_vector,
                                                bathymetry_gradient);

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
            double_contract(fe_values[momentum_extractor].gradient(i, q),
                            momentum_inviscid_flux[q] + momentum_viscous_flux[q])
            // bathymetry source term
            -
            fe_values[momentum_extractor].value(i, q) * gravity * height[q] *
              bathymetry_gradient[q]) *
          fe_values.JxW(q);
      }
    }

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
  std::vector<Tensor<1, dim>> velocity(n);
  compute_velocity(height, momentum, velocity);

  // loop over vector elements
  for (unsigned int q = 0; q < n; ++q)
  {
    // compute height inviscid flux
    height_flux[q] = momentum[q];

    // compute momentum inviscid flux
    Tensor<2, dim> velocity_times_momentum;
    outer_product(velocity_times_momentum, velocity[q], momentum[q]);
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
  std::vector<Tensor<1, dim>> & height_viscous_flux,
  std::vector<Tensor<2, dim>> & momentum_viscous_flux) const
{
  // get number of vector elements
  const unsigned int n = height_gradient.size();

  // loop over vector elements
  for (unsigned int q = 0; q < n; ++q)
  {
    // density viscous flux
    height_viscous_flux[q] = -viscosity * height_gradient[q];

    // momentum viscous flux
    momentum_viscous_flux[q] = -viscosity * momentum_gradient[q];
  }
}

/**
 * \brief Computes a vector of velocity values.
 *
 * \param[in] height  vector of height values
 * \param[in] momentum vector of momentum values
 * \param[out] velocity vector of velocity values
 */
template <int dim>
void ShallowWater<dim>::compute_velocity(
  const std::vector<double> & height,
  const std::vector<Tensor<1, dim>> & momentum,
  std::vector<Tensor<1, dim>> & velocity) const
{
  unsigned int n = height.size();
  for (unsigned int q = 0; q < n; ++q)
    velocity[q] = momentum[q] / height[q];
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
  cell_iterator cell = this->dof_handler.begin_active(),
                endc = this->dof_handler.end();
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
    std::vector<Tensor<1, dim>> velocity(this->n_q_points_cell);
    compute_velocity(height, momentum, velocity);

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
 * \brief Computes entropy \f$\eta\f$ at each degree of freedom point.
 *
 * For the shallow water equations, the entropy is defined as
 * \f[
 *   \eta(\mathbf{u}) = \frac{1}{2}\frac{\mathbf{q}\cdot\mathbf{q}}{h}
 *   + \frac{1}{2}g h^2
 * \f]
 */
/*
template <int dim>
void ShallowWater<dim>::compute_entropy_dof(const Vector<double> & solution,
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
*/

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
 * \brief Computes the divergence of the entropy flux
 *        \f$\nabla\cdot\mathbf{f}^\eta\f$ at each quadrature point
 *        on a cell or face.
 */
template <int dim>
void ShallowWater<dim>::compute_divergence_entropy_flux(
  const Vector<double> & solution,
  const FEValuesBase<dim> & fe_values,
  Vector<double> & divergence_entropy_flux) const
{
  // get number of quadrature points
  const unsigned int n = divergence_entropy_flux.size();

  // get solution values
  std::vector<double> height(n);
  std::vector<Tensor<1, dim>> momentum(n);
  fe_values[height_extractor].get_function_values(solution, height);
  fe_values[momentum_extractor].get_function_values(solution, momentum);

  // get height gradient
  std::vector<Tensor<1, dim>> height_gradient(n);
  fe_values[height_extractor].get_function_gradients(solution, height_gradient);

  // get momentum divergence
  std::vector<double> momentum_divergence(n);
  fe_values[momentum_extractor].get_function_divergences(solution,
                                                         momentum_divergence);

  // identity tensor
  SymmetricTensor<2, dim> identity_tensor_sym = unit_symmetric_tensor<dim>();
  Tensor<2, dim> identity_tensor(identity_tensor_sym);

  // compute divergence of entropy flux for each quadrature point
  for (unsigned int q = 0; q < n; ++q)
  {
    // compute derivative of entropy with respect to each component
    const double dentropy_dheight =
      -0.5 * momentum[q] * momentum[q] / (height[q] * height[q]) +
      gravity * height[q];
    const Tensor<1, dim> dentropy_dmomentum = momentum[q] / height[q];

    // compute outer product of momentum with itself
    Tensor<2, dim> momentum_times_momentum;
    outer_product(momentum_times_momentum, momentum[q], momentum[q]);

    // compute divergence of each component flux
    const Tensor<2, dim> aux = gravity * height[q] * identity_tensor;
    const double divergence_height_flux = momentum_divergence[q];
    const Tensor<1, dim> divergence_momentum_flux =
      2.0 * momentum[q] / height[q] * momentum_divergence[q] +
      (aux - momentum_times_momentum / (height[q] * height[q])) *
        height_gradient[q];

    // compute divergence of entropy flux
    divergence_entropy_flux[q] = dentropy_dheight * divergence_height_flux +
      dentropy_dmomentum * divergence_momentum_flux;
  }
}

template <int dim>
void ShallowWater<dim>::output_results(PostProcessor<dim> & postprocessor) const
{
  // create shallow water post-processor
  ShallowWaterPostProcessor<dim> shallowwater_postprocessor(bathymetry_function);

  // call post-processor
  postprocessor.output_results(this->new_solution,
                               this->dof_handler,
                               this->triangulation,
                               shallowwater_postprocessor);
}