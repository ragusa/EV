/**
 * \file Euler.cc
 * \brief Provides the function definitions for the Euler class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] params Euler equation parameters
 */
template <int dim>
Euler<dim>::Euler(const EulerParameters<dim> & params)
  : ConservationLaw<dim>(params),
    euler_parameters(params),
    density_extractor(0),
    momentum_extractor(1),
    energy_extractor(1 + dim)
{
}

template <int dim>
std::vector<std::string> Euler<dim>::get_component_names()
{
  std::vector<std::string> names(n_euler_components);

  names[0] = "density";
  for (int d = 0; d < dim; ++d)
    names[1 + d] = "momentum";
  names[dim + 1] = "energy";

  return names;
}

template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation> Euler<
  dim>::get_component_interpretations()
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretations(n_euler_components);

  component_interpretations[0] = DataComponentInterpretation::component_is_scalar;
  for (int d = 0; d < dim; ++d)
    component_interpretations[1 + d] =
      DataComponentInterpretation::component_is_part_of_vector;
  component_interpretations[dim + 1] =
    DataComponentInterpretation::component_is_scalar;

  return component_interpretations;
}

template <int dim>
void Euler<dim>::define_problem()
{
  switch (euler_parameters.problem_id)
  {
    case 0: // 1-D Sod shock tube problem
    {
      Assert(dim == 1, ExcImpossibleInDim(dim));

      // name of problem
      this->problem_name = "sod_tube";

      // create domain
      double domain_start = 0;
      double domain_width = 1.0;
      this->domain_volume = std::pow(domain_width, dim);
      GridGenerator::hyper_cube(
        this->triangulation, domain_start, domain_start + domain_width);

      // only 1 type of BC: zero Dirichlet
      this->n_boundaries = 1;

      // set all boundary indicators to zero
      typename Triangulation<dim>::cell_iterator cell =
                                                   this->triangulation.begin(),
                                                 endc = this->triangulation.end();
      for (; cell != endc; ++cell)
        for (unsigned int face = 0; face < this->faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_indicator(0);

      // set boundary conditions type for each boundary and component
      this->boundary_types.resize(this->n_boundaries);
      this->boundary_types[0].resize(this->n_components);
      this->boundary_types[0][0] =
        ConservationLaw<dim>::dirichlet; // density has Dirichlet BC
      this->boundary_types[0][1] =
        ConservationLaw<dim>::dirichlet; // x-momentum has Dirichlet BC
      this->boundary_types[0][2] =
        ConservationLaw<dim>::dirichlet; // energy has Dirichlet BC

      // set function strings to be parsed for dirichlet boundary condition
      // functions
      this->use_exact_solution_as_BC = false;
      this->dirichlet_function_strings.resize(this->n_boundaries);
      for (unsigned int boundary = 0; boundary < this->n_boundaries; ++boundary)
      {
        this->dirichlet_function_strings[boundary].resize(this->n_components);
        this->dirichlet_function_strings[boundary][0] =
          "if(x<0.5,1.0,0.125)";                             // BC for density
        this->dirichlet_function_strings[boundary][1] = "0"; // BC for x-momentum
        this->dirichlet_function_strings[boundary][2] =
          "if(x<0.5,2.5,0.25)"; // BC for energy
      }

      // initial conditions for each solution component
      this->initial_conditions_strings[0] = "if(x<0.5,1.0,0.125)";
      this->initial_conditions_strings[1] = "0";
      this->initial_conditions_strings[2] = "if(x<0.5,2.5,0.25)";

      // problem parameters
      const double rho_left = 1.0;
      const double u_left = 0.0;
      const double p_left = 1.0;
      const double rho_right = 0.125;
      const double u_right = 0.0;
      const double p_right = 0.1;
      gamma = 1.4;
      const double x_interface = 0.5;

      this->has_exact_solution = true;
      // create and initialize Riemann solver for exact solution
      std::shared_ptr<EulerRiemannSolver<dim>> exact_solution_function_derived =
        std::make_shared<EulerRiemannSolver<dim>>(rho_left,
                                                  u_left,
                                                  p_left,
                                                  rho_right,
                                                  u_right,
                                                  p_right,
                                                  gamma,
                                                  x_interface);
      // point base class pointer to derived class function object
      this->exact_solution_function = exact_solution_function_derived;

      break;
    }
    case 1: // 1-D Leblanc tube problem
    {
      Assert(dim == 1, ExcImpossibleInDim(dim));

      // name of problem
      this->problem_name = "leblanc_tube";

      // create domain
      double domain_start = 0;
      double domain_width = 1.0;
      this->domain_volume = std::pow(domain_width, dim);
      GridGenerator::hyper_cube(
        this->triangulation, domain_start, domain_start + domain_width);

      // only 1 type of BC: zero Dirichlet
      this->n_boundaries = 1;

      // set all boundary indicators to zero
      typename Triangulation<dim>::cell_iterator cell =
                                                   this->triangulation.begin(),
                                                 endc = this->triangulation.end();
      for (; cell != endc; ++cell)
        for (unsigned int face = 0; face < this->faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_indicator(0);

      // set boundary conditions type for each boundary and component
      this->boundary_types.resize(this->n_boundaries);
      this->boundary_types[0].resize(this->n_components);
      this->boundary_types[0][0] =
        ConservationLaw<dim>::dirichlet; // density has Dirichlet BC
      this->boundary_types[0][1] =
        ConservationLaw<dim>::dirichlet; // x-momentum has Dirichlet BC
      this->boundary_types[0][2] =
        ConservationLaw<dim>::dirichlet; // energy has Dirichlet BC

      // set function strings to be parsed for dirichlet boundary condition
      // functions
      this->use_exact_solution_as_BC = false;
      this->dirichlet_function_strings.resize(this->n_boundaries);
      for (unsigned int boundary = 0; boundary < this->n_boundaries; ++boundary)
      {
        this->dirichlet_function_strings[boundary].resize(this->n_components);
        this->dirichlet_function_strings[boundary][0] =
          "if(x<0.5,1.0,0.001)";                             // BC for density
        this->dirichlet_function_strings[boundary][1] = "0"; // BC for x-momentum
        this->dirichlet_function_strings[boundary][2] =
          "if(x<0.5,0.1,1.0e-10)"; // BC for energy
      }

      // initial conditions
      this->initial_conditions_strings[0] =
        "if(x<0.5,1.0,0.001)";                   // IC for density
      this->initial_conditions_strings[1] = "0"; // IC for x-momentum
      this->initial_conditions_strings[2] =
        "if(x<0.5,0.1,1.0e-10)"; // IC for energy

      // exact solution
      this->has_exact_solution = false;

      // physical constants
      gamma = 5.0 / 3.0;

      break;
    }
    case 2: // 2-D Noh Problem
    {
      // this is a 2-D problem
      Assert(dim == 2, ExcImpossibleInDim(dim));

      // name of problem
      this->problem_name = "noh";

      // create domain
      this->domain_volume =
        1.0; // domain is the unit hypercube, so domain volume is 1^dim
      GridIn<dim> input_grid;
      input_grid.attach_triangulation(this->triangulation);
      std::ifstream input_file("mesh/unit_square.msh");
      input_grid.read_msh(input_file);

      // four boundaries: each side of unit square
      this->n_boundaries = 4;

      // set boundary indicators
      double small_number = 1.0e-15;
      typename Triangulation<dim>::cell_iterator cell =
                                                   this->triangulation.begin(),
                                                 endc = this->triangulation.end();
      for (; cell != endc; ++cell)
        for (unsigned int face = 0; face < this->faces_per_cell; ++face)
          if (cell->face(face)->at_boundary())
          {
            Point<dim> face_center = cell->face(face)->center();
            if (face_center(1) < small_number)
            { // y = 0 boundary
              cell->face(face)->set_boundary_indicator(0);
            }
            else if (face_center(0) > 1.0 - small_number)
            { // x = 1 boundary
              cell->face(face)->set_boundary_indicator(1);
            }
            else if (face_center(1) > 1.0 - small_number)
            { // y = 1 boundary
              cell->face(face)->set_boundary_indicator(2);
            }
            else if (face_center(0) < small_number)
            { // x = 0 boundary
              cell->face(face)->set_boundary_indicator(3);
            }
            else
            {
              // all faces should have satisfied one of the conditions
              std::cout << "x = " << face_center(0) << std::endl;
              std::cout << "y = " << face_center(1) << std::endl;
              Assert(false, ExcInternalError());
            }
          }
      // set boundary conditions type for each boundary and component
      this->boundary_types.resize(this->n_boundaries);
      this->dirichlet_function_strings.resize(this->n_boundaries);
      for (unsigned int boundary = 0; boundary < this->n_boundaries; ++boundary)
      {
        this->boundary_types[boundary].resize(this->n_components);
        this->dirichlet_function_strings[boundary].resize(this->n_components);
        switch (boundary)
        {
          case 0: // y = 0 boundary
          {
            // all are reflective (zero Neumann)
            for (unsigned int component = 0; component < this->n_components;
                 ++component)
              this->boundary_types[boundary][component] =
                ConservationLaw<dim>::neumann;
            break;
          }
          case 1: // x = 1 boundary
          {
            // all are Dirichlet
            for (unsigned int component = 0; component < this->n_components;
                 ++component)
              this->boundary_types[boundary][component] =
                ConservationLaw<dim>::dirichlet;
            this->dirichlet_function_strings[boundary][0] =
              "if(sqrt(x^2+y^2)<t/3.0,16,1)"; // density
            this->dirichlet_function_strings[boundary][1] =
              "if(sqrt(x^2+y^2)<t/3.0,0,-x/sqrt(x^2+y^2))"; // mx
            this->dirichlet_function_strings[boundary][2] =
              "if(sqrt(x^2+y^2)<t/3.0,0,-y/sqrt(x^2+y^2))"; // my
            this->dirichlet_function_strings[boundary][3] =
              "if(sqrt(x^2+y^2)<t/3.0,16.0/3.0/(5.0/3.0-1),";
            this->dirichlet_function_strings[boundary][3] +=
              "1e-9/(5.0/3.0-1)+0.5)"; // energy
            break;
          }
          case 2: // y = 1 boundary
          {
            // all are Dirichlet
            for (unsigned int component = 0; component < this->n_components;
                 ++component)
              this->boundary_types[boundary][component] =
                ConservationLaw<dim>::dirichlet;
            this->dirichlet_function_strings[boundary][0] =
              "if(sqrt(x^2+y^2)<t/3.0,16,1)"; // density
            this->dirichlet_function_strings[boundary][1] =
              "if(sqrt(x^2+y^2)<t/3.0,0,-x/sqrt(x^2+y^2))"; // mx
            this->dirichlet_function_strings[boundary][2] =
              "if(sqrt(x^2+y^2)<t/3.0,0,-y/sqrt(x^2+y^2))"; // my
            this->dirichlet_function_strings[boundary][3] =
              "if(sqrt(x^2+y^2)<t/3.0,16.0/3.0/(5.0/3.0-1),";
            this->dirichlet_function_strings[boundary][3] +=
              "1e-9/(5.0/3.0-1)+0.5)"; // energy
            break;
          }
          case 3: // x = 0 boundary
          {
            // all are reflective (zero Neumann)
            for (unsigned int component = 0; component < this->n_components;
                 ++component)
              this->boundary_types[boundary][component] =
                ConservationLaw<dim>::neumann;
            break;
          }
          default:
          {
            std::cout << "boundary indicator is " << boundary << std::endl;
            Assert(false, ExcInternalError());
            break;
          }
        }
      }
      this->use_exact_solution_as_BC = false;

      // initial conditions for each solution component
      this->initial_conditions_strings[0] = "1";
      this->initial_conditions_strings[1] = "if(x==0,0,-x/sqrt(x^2+y^2))";
      this->initial_conditions_strings[2] = "if(y==0,0,-y/sqrt(x^2+y^2))";
      this->initial_conditions_strings[3] = "1e-9/(5.0/3.0-1)+0.5";

      // exact solution
      this->has_exact_solution = true;
      this->exact_solution_strings[0] = "if(sqrt(x^2+y^2)<t/3.0,16,1)"; // density
      this->exact_solution_strings[1] =
        "if(sqrt(x^2+y^2)<t/3.0,0,-x/sqrt(x^2+y^2))"; // mx
      this->exact_solution_strings[2] =
        "if(sqrt(x^2+y^2)<t/3.0,0,-y/sqrt(x^2+y^2))"; // my
      this->exact_solution_strings[3] =
        "if(sqrt(x^2+y^2)<t/3.0,16.0/3.0/(5.0/3.0-1),";
      this->exact_solution_strings[3] += "1e-9/(5.0/3.0-1)+0.5)"; // energy

      // physical constants
      gamma = 5.0 / 3.0;

      break;
    }
    default:
    {
      Assert(false, ExcNotImplemented());
      break;
    }
  }
}

template <int dim>
void Euler<dim>::assemble_lumped_mass_matrix()
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
          local_mass(i, i) += (fe_values[density_extractor].value(i, q) *
                                 fe_values[density_extractor].value(j, q) +
                               fe_values[momentum_extractor].value(i, q) *
                                 fe_values[momentum_extractor].value(j, q) +
                               fe_values[energy_extractor].value(i, q) *
                                 fe_values[energy_extractor].value(j, q)) *
            fe_values.JxW(q);
        }

    // add to global mass matrix with contraints
    this->constraints.distribute_local_to_global(
      local_mass, local_dof_indices, this->lumped_mass_matrix);
  }
}

template <int dim>
void Euler<dim>::compute_ss_residual(Vector<double> & f)
{
  // reset vector
  f = 0.0;

  FEValues<dim> fe_values(this->fe,
                          this->cell_quadrature,
                          update_values | update_gradients | update_JxW_values);
  FEFaceValues<dim> fe_face_values(this->fe,
                                   this->face_quadrature,
                                   update_values | update_normal_vectors |
                                     update_JxW_values);

  Vector<double> cell_residual(this->dofs_per_cell);
  std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);

  // loop over cells
  cell_iterator cell = this->dof_handler.begin_active(),
                endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // reset cell residual
    cell_residual = 0;

    // compute local residual
    compute_cell_ss_residual(fe_values, fe_face_values, cell, cell_residual);

    // aggregate local residual into global residual
    cell->get_dof_indices(local_dof_indices);
    this->constraints.distribute_local_to_global(
      cell_residual, local_dof_indices, f);
  } // end cell loop
}

/** \brief Computes the contribution of the steady-state residual
 *         from the current cell, including face contributions.
 *  \param[in] fe_values FEValues object
 *  \param[in] fe_face_values FEFaceValues object
 *  \param[in] cell cell iterator for current cell
 *  \param[out] cell_residual cell contribution to global right hand side
 */
template <int dim>
void Euler<dim>::compute_cell_ss_residual(FEValues<dim> & fe_values,
                                          FEFaceValues<dim> & fe_face_values,
                                          const cell_iterator & cell,
                                          Vector<double> & cell_residual)
{
  // reinitialize fe values for cell
  fe_values.reinit(cell);

  // get solution values
  std::vector<double> density(this->n_q_points_cell);
  std::vector<Tensor<1, dim>> momentum(this->n_q_points_cell);
  std::vector<double> energy(this->n_q_points_cell);
  fe_values[density_extractor].get_function_values(this->new_solution, density);
  fe_values[momentum_extractor].get_function_values(this->new_solution, momentum);
  fe_values[energy_extractor].get_function_values(this->new_solution, energy);

  // get solution gradients
  std::vector<Tensor<1, dim>> density_gradient(this->n_q_points_cell);
  std::vector<Tensor<2, dim>> momentum_gradient(this->n_q_points_cell);
  std::vector<Tensor<1, dim>> energy_gradient(this->n_q_points_cell);
  fe_values[density_extractor].get_function_gradients(this->new_solution,
                                                      density_gradient);
  fe_values[momentum_extractor].get_function_gradients(this->new_solution,
                                                       momentum_gradient);
  fe_values[energy_extractor].get_function_gradients(this->new_solution,
                                                     energy_gradient);

  // get momentum divergences
  std::vector<double> momentum_divergence(this->n_q_points_cell);
  fe_values[momentum_extractor].get_function_divergences(this->new_solution,
                                                         momentum_divergence);

  // get viscosity values on cell
  std::vector<double> viscosity(this->n_q_points_cell);
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    viscosity[q] = this->viscosity[cell];

  // compute inviscid fluxes
  std::vector<Tensor<1, dim>> density_inviscid_flux(this->n_q_points_cell);
  std::vector<Tensor<2, dim>> momentum_inviscid_flux(this->n_q_points_cell);
  std::vector<Tensor<1, dim>> energy_inviscid_flux(this->n_q_points_cell);
  compute_inviscid_fluxes(density,
                          momentum,
                          energy,
                          density_inviscid_flux,
                          momentum_inviscid_flux,
                          energy_inviscid_flux);

  // compute viscous fluxes
  std::vector<Tensor<1, dim>> density_viscous_flux(this->n_q_points_cell);
  std::vector<Tensor<2, dim>> momentum_viscous_flux(this->n_q_points_cell);
  std::vector<Tensor<1, dim>> energy_viscous_flux(this->n_q_points_cell);
  compute_viscous_fluxes(viscosity,
                         density,
                         momentum,
                         energy,
                         density_gradient,
                         momentum_gradient,
                         energy_gradient,
                         momentum_divergence,
                         density_viscous_flux,
                         momentum_viscous_flux,
                         energy_viscous_flux);

  // loop over quadrature points
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
  {
    // loop over test functions
    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      // sum up terms into cell residual
      cell_residual(i) +=
        ((
           // density
           fe_values[density_extractor].gradient(i, q) *
             (density_inviscid_flux[q] + density_viscous_flux[q])
           // momentum
           +
           double_contract(fe_values[momentum_extractor].gradient(i, q),
                           momentum_inviscid_flux[q] + momentum_viscous_flux[q])
           // energy
           +
           fe_values[energy_extractor].gradient(i, q) *
             (energy_inviscid_flux[q] + energy_viscous_flux[q])) *
         fe_values.JxW(q));
  }

  // resize flux vectors for face computations
  density_inviscid_flux.resize(this->n_q_points_face);
  momentum_inviscid_flux.resize(this->n_q_points_face);
  energy_inviscid_flux.resize(this->n_q_points_face);

  // loop over faces
  for (unsigned int face = 0; face < this->faces_per_cell; ++face)
  {
    // add term for boundary faces
    if (cell->at_boundary(face))
    {
      fe_face_values.reinit(cell, face);

      // get values on face
      std::vector<double> density(this->n_q_points_face);
      std::vector<Tensor<1, dim>> momentum(this->n_q_points_face);
      std::vector<double> energy(this->n_q_points_face);
      fe_face_values[density_extractor].get_function_values(this->new_solution,
                                                            density);
      fe_face_values[momentum_extractor].get_function_values(this->new_solution,
                                                             momentum);
      fe_face_values[energy_extractor].get_function_values(this->new_solution,
                                                           energy);

      // compute inviscid fluxes
      compute_inviscid_fluxes(density,
                              momentum,
                              energy,
                              density_inviscid_flux,
                              momentum_inviscid_flux,
                              energy_inviscid_flux);

      /*
               std::vector<Tensor<1, dim> >
         solution_gradients_face(this->n_q_points_face);
               fe_face_values[velocity].get_function_gradients
         (this->new_solution,solution_gradients_face);

               // compute viscosity
               Vector<double> viscosity_face(this->n_q_points_face);
               switch (euler_parameters.viscosity_type)
               {
                  case EulerParameters<dim>::constant:
                  {
                     for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                        viscosity_face(q) =
         euler_parameters.constant_viscosity_value;
                     break;
                  }
                  case EulerParameters<dim>::first_order:
                  {
                     // get max velocity on cell
                     std::vector<double> local_solution(this->n_q_points_face);
                     fe_face_values.get_function_values(this->new_solution,
         local_solution);
                     double max_velocity = 0.0;
                     for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                        max_velocity = std::max( max_velocity, local_solution[q]);

                     // compute first-order viscosity
                     double cell_diameter = cell->diameter();
                     double viscosity_value = 0.5 *
         euler_parameters.first_order_viscosity_coef * cell_diameter *
         max_velocity;
                     for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                        viscosity_face(q) = viscosity_value;

                     break;
                  }
                  default:
                  {
                     Assert(false,ExcNotImplemented());
                     break;
                  }
               }
      */
      // loop over test functions
      for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
        // loop over quadrature points
        for (unsigned int q = 0; q < this->n_q_points_face; ++q)
          cell_residual(i) -= (
                                // density
                                fe_face_values[density_extractor].value(i, q) *
                                  density_inviscid_flux[q]
                                // momentum
                                +
                                fe_face_values[momentum_extractor].value(i, q) *
                                  momentum_inviscid_flux[q]
                                // energy
                                +
                                fe_face_values[energy_extractor].value(i, q) *
                                  energy_inviscid_flux[q]
                                // normal vector and JxW
                                ) *
            fe_face_values.normal_vector(q) * fe_face_values.JxW(q);
    }
  }
}

/**
 * \brief Computes the viscous fluxes.
 */
template <int dim>
void Euler<dim>::compute_viscous_fluxes(
  const std::vector<double> & viscosity,
  const std::vector<double> & density,
  const std::vector<Tensor<1, dim>> & momentum,
  const std::vector<double> & energy,
  const std::vector<Tensor<1, dim>> & density_gradient,
  const std::vector<Tensor<2, dim>> & momentum_gradient,
  const std::vector<Tensor<1, dim>> & energy_gradient,
  const std::vector<double> & momentum_divergence,
  std::vector<Tensor<1, dim>> & density_viscous_flux,
  std::vector<Tensor<2, dim>> & momentum_viscous_flux,
  std::vector<Tensor<1, dim>> & energy_viscous_flux) const
{
  // get number of vector elements
  const unsigned int n = density.size();

  // compute auxiliary quantities
  std::vector<Tensor<1, dim>> velocity(n);
  std::vector<double> pressure(n);
  std::vector<double> internal_energy(n);
  std::vector<double> temperature(n);
  compute_velocity(density, momentum, velocity);
  compute_pressure(density, momentum, energy, pressure);
  compute_internal_energy(density, momentum, energy, internal_energy);
  compute_temperature(internal_energy, temperature);

  // loop over vector elements
  for (unsigned int q = 0; q < n; ++q)
  {
    // density viscous flux
    double nu = viscosity[q];
    density_viscous_flux[q] = -nu * density_gradient[q];

    // momentum viscous flux
    Tensor<2, dim> density_viscous_flux_times_velocity;
    outer_product(
      density_viscous_flux_times_velocity, density_viscous_flux[q], velocity[q]);
    Tensor<2, dim> momentum_times_density_gradient;
    outer_product(
      momentum_times_density_gradient, momentum[q], density_gradient[q]);
    Tensor<2, dim> velocity_symmetric_gradient =
      (momentum_gradient[q] * density[q] - momentum_times_density_gradient) /
      (density[q] * density[q]);
    momentum_viscous_flux[q] = -nu * density[q] * velocity_symmetric_gradient +
      density_viscous_flux_times_velocity;

    // energy viscous flux
    double kappa = nu;
    double sp = -pressure[q] / (temperature[q] * density[q] * density[q]);
    double se = 1.0 / temperature[q];
    double velocity_divergence =
      (momentum_divergence[q] * density[q] - momentum[q] * density_gradient[q]) /
      (density[q] * density[q]);
    Tensor<1, dim> pressure_gradient =
      (gamma - 1.0) * (energy_gradient[q] - velocity_divergence * momentum[q] -
                       0.5 * velocity[q] * velocity[q] * density_gradient[q]);
    Tensor<1, dim> internal_energy_gradient =
      (pressure_gradient * density[q] - pressure[q] * density_gradient[q]) /
      ((gamma - 1.0) * density[q] * density[q]);
    Tensor<1, dim> l_flux =
      (nu - kappa) * density[q] * sp / se * density_gradient[q] -
      nu * internal_energy[q] * density_gradient[q] -
      kappa * density[q] * internal_energy_gradient;
    Tensor<1, dim> h_flux =
      l_flux - 0.5 * velocity[q] * velocity[q] * density_viscous_flux[q];
    energy_viscous_flux[q] = h_flux + momentum_viscous_flux[q] * velocity[q];
  }
}

/** \brief Computes the contribution of the steady-state residual
 *         from the faces of the current cell.
 *  \param [in] fe_face_values FEFaceValues object
 *  \param [in] cell cell iterator for current cell
 *  \param [out] cell_residual cell contribution to global right hand side
 */
/*
template <int dim>
void Euler<dim>::compute_face_ss_residual(FEFaceValues<dim> &fe_face_values,
                                          const cell_iterator &cell,
                                          Vector<double> &cell_residual)
{
   // loop over faces
   for (unsigned int face = 0; face < this->faces_per_cell; ++face)
   {
      // add term for boundary faces
      if (cell->at_boundary(face))
      {
         fe_face_values.reinit(cell, face);

         std::vector<Tensor<1, dim> >
solution_gradients_face(this->n_q_points_face);
         fe_face_values[velocity].get_function_gradients
(this->new_solution,solution_gradients_face);

         // compute viscosity
         Vector<double> viscosity_face(this->n_q_points_face);
         switch (euler_parameters.viscosity_type)
         {
            case EulerParameters<dim>::constant:
            {
               for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                  viscosity_face(q) = euler_parameters.constant_viscosity_value;
               break;
            }
            case EulerParameters<dim>::first_order:
            {
               // get max velocity on cell
               std::vector<double> local_solution(this->n_q_points_face);
               fe_face_values.get_function_values(this->new_solution,
local_solution);
               double max_velocity = 0.0;
               for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                  max_velocity = std::max( max_velocity, local_solution[q]);

               // compute first-order viscosity
               double cell_diameter = cell->diameter();
               double viscosity_value = 0.5 *
euler_parameters.first_order_viscosity_coef * cell_diameter * max_velocity;
               for (unsigned int q = 0; q < this->n_q_points_face; ++q)
                  viscosity_face(q) = viscosity_value;

               break;
            }
            default:
            {
               Assert(false,ExcNotImplemented());
               break;
            }
         }
         // loop over test functions
         for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
            // loop over quadrature points
            for (unsigned int q = 0; q < this->n_q_points_face; ++q)
               cell_residual(i) += fe_face_values[velocity].value(i,q)
    *viscosity_face(q)
    *solution_gradients_face[q]
    *fe_face_values.normal_vector(q)
    *fe_face_values.JxW(q);
      }
   }
}
*/

/**
 * \brief Computes the inviscid fluxes required to be evaluated in cell and face
 *        integrations.
 *
 * \param[in] density  vector of density values
 * \param[in] momentum vector of momentum values
 * \param[in] energy   vector of energy values
 * \param[out] density_flux  vector of density inviscid flux values
 * \param[out] momentum_flux vector of momentum inviscid flux values
 * \param[out] energy_flux   vector of energy inviscid flux values
*/
template <int dim>
void Euler<dim>::compute_inviscid_fluxes(
  const std::vector<double> & density,
  const std::vector<Tensor<1, dim>> & momentum,
  const std::vector<double> & energy,
  std::vector<Tensor<1, dim>> & density_flux,
  std::vector<Tensor<2, dim>> & momentum_flux,
  std::vector<Tensor<1, dim>> & energy_flux) const
{
  // get number of vector elements
  const unsigned int n = density.size();

  // identity tensor
  SymmetricTensor<2, dim> identity_tensor = unit_symmetric_tensor<dim>();

  // compute auxiliary quantities
  std::vector<Tensor<1, dim>> velocity(n);
  std::vector<double> pressure(n);
  compute_velocity(density, momentum, velocity);
  compute_pressure(density, momentum, energy, pressure);

  // loop over vector elements
  for (unsigned int q = 0; q < n; ++q)
  {
    // compute density inviscid flux
    density_flux[q] = momentum[q];

    // compute momentum inviscid flux
    Tensor<2, dim> velocity_times_momentum;
    outer_product(velocity_times_momentum, velocity[q], momentum[q]);
    momentum_flux[q] = velocity_times_momentum + pressure[q] * identity_tensor;

    // compute energy inviscid flux
    energy_flux[q] = velocity[q] * (energy[q] + pressure[q]);
  }
}

/** \brief computes the steady-state Jacobian, which is used if an implicit
 *         time integration scheme is to be used.
 */
/*
template <int dim>
void Euler<dim>::compute_ss_jacobian()
{
   // only been coded in 1-D
   Assert(dim==1,ExcNotImplemented());

   // get current solution values and gradients
   std::vector<double>         density          (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > density_gradient (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > momentum         (this->n_q_points_cell);
   std::vector<Tensor<2,dim> > momentum_gradient(this->n_q_points_cell);
   std::vector<SymmetricTensor<2,dim> >
momentum_symmetric_gradient(this->n_q_points_cell);
   std::vector<double>         momentum_divergence (this->n_q_points_cell);
   std::vector<double>         energy           (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > energy_gradient  (this->n_q_points_cell);
   std::vector<double>         internal_energy  (this->n_q_points_cell);
   std::vector<double>         temperature      (this->n_q_points_cell);
   std::vector<Tensor<1,dim> > velocity         (this->n_q_points_cell);
   std::vector<double>         pressure         (this->n_q_points_cell);
   std::vector<double>         dpdrho           (this->n_q_points_cell);
   std::vector<double>         dpdmx            (this->n_q_points_cell);
   std::vector<double>         dpdE             (this->n_q_points_cell);

   // derivatives of Euler flux functions
   Tensor<1,dim> dfdu_rho_rho;
   Tensor<1,dim> dfdu_rho_mx;
   Tensor<1,dim> dfdu_rho_E;
   Tensor<1,dim> dfdu_mx_rho;
   Tensor<1,dim> dfdu_mx_mx;
   Tensor<1,dim> dfdu_mx_E;
   Tensor<1,dim> dfdu_E_rho;
   Tensor<1,dim> dfdu_E_mx;
   Tensor<1,dim> dfdu_E_E;

   // other
   Tensor<1,dim> unit_vector_x; unit_vector_x[0] = 1.0;

   // local DoF indices
   std::vector<unsigned int> local_dof_indices(this->dofs_per_cell);

   // cell matrix
   FullMatrix<double> cell_matrix(this->dofs_per_cell,this->dofs_per_cell);

   // FE values
   FEValues<dim> fe_values(this->fe,this->cell_quadrature,
      update_values | update_gradients | update_JxW_values);

   // ones vector
   Tensor<1,dim> ones_vector;
   for (unsigned int d = 0; d < dim; ++d)
      ones_vector[d] = 1.0;

   // reset steady-state Jacobian to zero
   this->system_matrix = 0.0;
   // loop over cells
   cell_iterator cell = this->dof_handler.begin_active(),
      endc = this->dof_handler.end();
   for (; cell!=endc; ++cell)
   {
      fe_values.reinit(cell);

      fe_values[density_extractor] .get_function_values   (this->new_solution,
density);
      fe_values[density_extractor] .get_function_gradients(this->new_solution,
density_gradient);
      fe_values[momentum_extractor].get_function_values   (this->new_solution,
momentum);
      fe_values[momentum_extractor].get_function_gradients(this->new_solution,
momentum_gradient);
      fe_values[momentum_extractor].get_function_symmetric_gradients
(this->new_solution, momentum_symmetric_gradient);
      fe_values[momentum_extractor].get_function_divergences(this->new_solution,
momentum_divergence);
      fe_values[energy_extractor]  .get_function_values     (this->new_solution,
energy);
      fe_values[energy_extractor]  .get_function_gradients  (this->new_solution,
energy_gradient);
      compute_pressure_cell (pressure, dpdrho, dpdmx, dpdE, density, momentum,
energy);
      dfdu_rho_rho = 0.0;
      dfdu_rho_mx  = unit_vector_x;
      dfdu_rho_E   = 0.0;
      dfdu_mx_rho  = 0.0;

      // reset cell matrix to zero
      cell_matrix = 0;
      // loop over quadrature points in cell
      for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
      {
         double rho = density[q];
         double mx = momentum[q]*unit_vector_x;
         double E = energy[q];
         double p = pressure[q];

         dfdu_mx_mx[0] = 2.0*mx/rho + dpdmx[q];
         dfdu_mx_E[0]  = dpdE[q];
         dfdu_E_rho[0] = -mx/std::pow(rho,2)*(E+p) + mx/rho*dpdrho[q];
         dfdu_E_mx[0]  = (E+p)/rho + mx/rho*dpdmx[q];
         dfdu_E_E[0]   = dpdE[q];

         for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
         {
            for (unsigned int j = 0; j < this->dofs_per_cell; ++j)
            {
               cell_matrix(i,j) += (
                  // conservation of mass
                  fe_values[density_extractor].gradient(i,q) * dfdu_rho_rho *
                  fe_values[density_extractor].value(j,q)
                  + fe_values[density_extractor].gradient(i,q) * dfdu_rho_mx *
                  //fe_values[momentum_extractor].value(j,q)
                  fe_values.shape_value_component(j,q,1)
                  + fe_values[density_extractor].gradient(i,q) * dfdu_rho_E *
                  fe_values[energy_extractor].value(j,q)

                  // conservation of x-momentum
                  //+ fe_values[momentum_extractor].shape_grad_component(i,q,0) *
dfdu_mx_rho *
                  + fe_values.shape_grad_component(i,q,1) * dfdu_mx_rho *
                  fe_values[density_extractor].value(j,q)
                  //+ fe_values[momentum_extractor].shape_grad_component(i,q,0) *
dfdu_mx_mx *
                  + fe_values.shape_grad_component(i,q,1) * dfdu_mx_mx *
                  fe_values.shape_value_component(j,q,1)
                  //+ fe_values[momentum_extractor].shape_grad_component(i,q,0) *
dfdu_mx_E *
                  + fe_values.shape_grad_component(i,q,1) * dfdu_mx_E *
                  fe_values[energy_extractor].value(j,q)

                  // conservation of energy
                  + fe_values[energy_extractor].gradient(i,q) * dfdu_E_rho *
                  fe_values[density_extractor].value(j,q)
                  + fe_values[energy_extractor].gradient(i,q) * dfdu_E_mx *
                  //fe_values[momentum_extractor].value(j,q)
                  fe_values.shape_value_component(j,q,1)
                  + fe_values[energy_extractor].gradient(i,q) * dfdu_E_E *
                  fe_values[energy_extractor].value(j,q)

               ) * fe_values.JxW(q);
            }
         }
      }
      // get dof indices
      cell->get_dof_indices(local_dof_indices);
      // aggregate cell matrix into global matrix
      this->constraints.distribute_local_to_global(cell_matrix, local_dof_indices,
this->system_matrix);
   }
}
*/

/**
 * \brief Computes a vector of velocity values.
 *
 * \param[in] density  vector of density values
 * \param[in] momentum vector of momentum values
 * \param[out] velocity vector of velocity values
 */
template <int dim>
void Euler<dim>::compute_velocity(const std::vector<double> & density,
                                  const std::vector<Tensor<1, dim>> & momentum,
                                  std::vector<Tensor<1, dim>> & velocity) const
{
  unsigned int n = density.size();
  for (unsigned int q = 0; q < n; ++q)
    velocity[q] = momentum[q] / density[q];
}

/**
 * \brief Computes a vector of internal energy values.
 */
template <int dim>
void Euler<dim>::compute_internal_energy(
  const std::vector<double> & density,
  const std::vector<Tensor<1, dim>> & momentum,
  const std::vector<double> & energy,
  std::vector<double> & internal_energy) const
{
  const unsigned int n = density.size();

  for (unsigned int q = 0; q < n; ++q)
    internal_energy[q] =
      (energy[q] - 0.5 * momentum[q] * momentum[q] / density[q]) / density[q];
}

/** \brief Computes internal energy at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_internal_energy_cell(
  std::vector<double> & internal_energy,
  const std::vector<double> & density,
  const std::vector<Tensor<1, dim>> & momentum,
  const std::vector<double> & energy) const
{
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    internal_energy[q] =
      (energy[q] - 0.5 * momentum[q] * momentum[q] / density[q]) / density[q];
}

/** \brief Computes internal energy at each quadrature point in face.
 */
template <int dim>
void Euler<dim>::compute_internal_energy_face(
  std::vector<double> & internal_energy,
  const std::vector<double> & density,
  const std::vector<Tensor<1, dim>> & momentum,
  const std::vector<double> & energy) const
{
  for (unsigned int q = 0; q < this->n_q_points_face; ++q)
    internal_energy[q] =
      (energy[q] - 0.5 * momentum[q] * momentum[q] / density[q]) / density[q];
}

/**
 * \brief Computes a vector of temperature values.
 *
 * \param[in] internal_energy internal energy
 * \param[in] temperature temperature
 */
template <int dim>
void Euler<dim>::compute_temperature(const std::vector<double> & internal_energy,
                                     std::vector<double> & temperature) const
{
  const unsigned int n = internal_energy.size();

  for (unsigned int q = 0; q < n; ++q)
    temperature[q] = (gamma - 1.0) * internal_energy[q];
}

/** \brief Computes temperature at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_temperature_cell(
  std::vector<double> & temperature,
  const std::vector<double> & internal_energy) const
{
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    temperature[q] = (gamma - 1.0) * internal_energy[q];
}

/** \brief Computes temperature at each quadrature point in face.
 */
template <int dim>
void Euler<dim>::compute_temperature_face(
  std::vector<double> & temperature,
  const std::vector<double> & internal_energy) const
{
  for (unsigned int q = 0; q < this->n_q_points_face; ++q)
    temperature[q] = (gamma - 1.0) * internal_energy[q];
}

/**
 * \brief Computes a vector of pressure values.
 *
 * \param[in] density  density
 * \param[in] momentum momentum
 * \param[in] energy   energy
 * \param[out] pressure pressure
 */
template <int dim>
void Euler<dim>::compute_pressure(const std::vector<double> & density,
                                  const std::vector<Tensor<1, dim>> & momentum,
                                  const std::vector<double> & energy,
                                  std::vector<double> & pressure) const
{
  unsigned int n = density.size();

  for (unsigned int q = 0; q < n; ++q)
    pressure[q] =
      (gamma - 1.0) * (energy[q] - 0.5 * momentum[q] * momentum[q] / density[q]);
}

/** \brief Computes pressure and pressure derivatives at each quadrature point in
 * cell.
 *  \param[out] pressure pressure
 *  \param[out] dpdrho derivative of pressure with respect to density
 *  \param[out] dpdmx derivative of pressure with respect to x-momentum
 *  \param[out] dpdE derivative of pressure with respect to energy
 *  \param[in] density density
 *  \param[in] momentum momentum
 *  \param[in] energy energy
 */
template <int dim>
void Euler<dim>::compute_pressure_cell(
  std::vector<double> & pressure,
  std::vector<double> & dpdrho,
  std::vector<double> & dpdmx,
  std::vector<double> & dpdE,
  const std::vector<double> & density,
  const std::vector<Tensor<1, dim>> & momentum,
  const std::vector<double> & energy) const
{
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
  {
    pressure[q] =
      (gamma - 1.0) * (energy[q] - 0.5 * momentum[q] * momentum[q] / density[q]);
    dpdrho[q] = 0.0;
    dpdmx[q] = (1.0 - gamma) * momentum[q][0];
    dpdE[q] = gamma - 1.0;
  }
}

/**
 * \brief Computes pressure and pressure derivatives at each quadrature point in
 * face.
 *
 * \param[out] pressure pressure
 * \param[out] dpdrho derivative of pressure with respect to density
 * \param[out] dpdmx derivative of pressure with respect to x-momentum
 * \param[out] dpdE derivative of pressure with respect to energy
 * \param[in] density density
 * \param[in] momentum momentum
 * \param[in] energy energy
 */
template <int dim>
void Euler<dim>::compute_pressure_face(
  std::vector<double> & pressure,
  std::vector<double> & dpdrho,
  std::vector<double> & dpdmx,
  std::vector<double> & dpdE,
  const std::vector<double> & density,
  const std::vector<Tensor<1, dim>> & momentum,
  const std::vector<double> & energy) const
{
  for (unsigned int q = 0; q < this->n_q_points_face; ++q)
  {
    pressure[q] = (gamma - 1.0) * (energy[q] - 0.5 * momentum[q] * momentum[q]);
    dpdrho[q] = 0.0;
    dpdmx[q] = (1.0 - gamma) * momentum[q][0];
    dpdE[q] = gamma - 1.0;
  }
}

/** \brief Computes speed of sound at each quadrature point in cell.
 */
template <int dim>
void Euler<dim>::compute_speed_of_sound(
  std::vector<double> & speed_of_sound,
  const std::vector<double> & density,
  const std::vector<double> & pressure) const
{
  for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    speed_of_sound[q] = std::sqrt(gamma * pressure[q] / density[q]);
}

template <int dim>
void Euler<dim>::update_flux_speeds()
{
  FEValues<dim> fe_values(this->fe, this->cell_quadrature, update_values);
  std::vector<double> density(this->n_q_points_cell);
  std::vector<Tensor<1, dim>> momentum(this->n_q_points_cell);
  std::vector<double> energy(this->n_q_points_cell);
  std::vector<Tensor<1, dim>> velocity(this->n_q_points_cell);
  std::vector<double> speed_of_sound(this->n_q_points_cell);
  std::vector<double> pressure(this->n_q_points_cell);

  // reset max flux speed
  this->max_flux_speed = 0.0;

  // loop over cells to compute first order viscosity at each quadrature point
  cell_iterator cell = this->dof_handler.begin_active(),
                endc = this->dof_handler.end();
  for (; cell != endc; ++cell)
  {
    // get conservative variables
    fe_values.reinit(cell);
    fe_values[density_extractor].get_function_values(this->new_solution, density);
    fe_values[momentum_extractor].get_function_values(this->new_solution,
                                                      momentum);
    fe_values[energy_extractor].get_function_values(this->new_solution, energy);

    // compute velocity
    compute_velocity(density, momentum, velocity);

    // compute speed of sound
    compute_pressure(density, momentum, energy, pressure);
    compute_speed_of_sound(speed_of_sound, density, pressure);

    // get maximum speed of sound and maximum fluid speed
    double max_speed_of_sound = 0.0;
    double max_speed = 0.0;
    for (unsigned int q = 0; q < this->n_q_points_cell; ++q)
    {
      max_speed_of_sound = std::max(max_speed_of_sound, speed_of_sound[q]);
      max_speed = std::max(max_speed, velocity[q].norm());
    }

    // get max flux speed
    this->max_flux_speed_cell[cell] = max_speed_of_sound + max_speed;
    this->max_flux_speed =
      std::max(this->max_flux_speed, this->max_flux_speed_cell[cell]);
  }
}

template <int dim>
void Euler<dim>::compute_entropy(const Vector<double> & solution,
                                 const FEValuesBase<dim> & fe_values,
                                 Vector<double> & entropy) const
{
  // get number of quadrature points
  const unsigned int n = entropy.size();

  std::vector<double> density(n);
  std::vector<Tensor<1, dim>> momentum(n);
  std::vector<double> energy(n);
  std::vector<double> pressure(n);

  // get conservative variables
  fe_values[density_extractor].get_function_values(this->new_solution, density);
  fe_values[momentum_extractor].get_function_values(this->new_solution, momentum);
  fe_values[energy_extractor].get_function_values(this->new_solution, energy);

  // compute pressure
  compute_pressure(density, momentum, energy, pressure);

  // compute entropy
  for (unsigned int q = 0; q < n; ++q)
    entropy[q] = density[q] / (gamma - 1.0) *
      std::log(pressure[q] / std::pow(density[q], gamma));
}

template <int dim>
void Euler<dim>::compute_divergence_entropy_flux(
  const Vector<double> & solution,
  const FEValuesBase<dim> & fe_values,
  Vector<double> & divergence) const
{
  // get number of quadrature points
  const unsigned int n = divergence.size();

  std::vector<double> density(n);
  std::vector<Tensor<1, dim>> momentum(n);
  std::vector<double> energy(n);
  std::vector<double> pressure(n);
  std::vector<Tensor<1, dim>> density_gradient(n);
  std::vector<Tensor<2, dim>> momentum_gradient(n);
  std::vector<Tensor<1, dim>> energy_gradient(n);
  std::vector<double> momentum_divergence(n);

  // get conservative variables
  fe_values[density_extractor].get_function_values(solution, density);
  fe_values[momentum_extractor].get_function_values(solution, momentum);
  fe_values[energy_extractor].get_function_values(solution, energy);

  // get gradients of conservative variables
  fe_values[density_extractor].get_function_gradients(solution,
                                                      density_gradient);
  fe_values[momentum_extractor].get_function_gradients(solution,
                                                       momentum_gradient);
  fe_values[energy_extractor].get_function_gradients(solution,
                                                     energy_gradient);

  // get divergence of momentum
  fe_values[momentum_extractor].get_function_divergences(solution,
                                                         momentum_divergence);

  // compute pressure
  compute_pressure(density, momentum, energy, pressure);

  // compute divergence of entropy flux
  for (unsigned int q = 0; q < n; ++q)
  {
    Tensor<1, dim> pressure_gradient =
      (gamma - 1.0) * energy_gradient[q] - momentum[q] * momentum_gradient[q];
    divergence[q] = momentum_divergence[q] *
        std::log(pressure[q] / std::pow(density[q], gamma)) +
      momentum[q] * (pressure_gradient / pressure[q] -
                     gamma / density[q] * density_gradient[q]);
  }
}
