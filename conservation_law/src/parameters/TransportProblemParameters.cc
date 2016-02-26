/**
 * \file TransportProblemParameters.cc
 * \brief Provides the function definitions for the TransportProblemParameters
 *        class.
 */

/**
 * \brief Constructor.
 *
 * \param[in] problem_name_  name of problem
 */
template <int dim>
TransportProblemParameters<dim>::TransportProblemParameters(
  const std::string & problem_name_)
  : ProblemParameters<dim>(problem_name_),
    cross_section_function(1),
    source_function(1)
{
  this->n_components = 1;
}

/**
 * \brief Declares derived parameters.
 */
template <int dim>
void TransportProblemParameters<dim>::declare_derived_parameters()
{
  // physics
  this->parameter_handler.enter_subsection("physics");
  {
    this->parameter_handler.declare_entry(
      "transport speed", "1.0", Patterns::Double(), "Transport speed");
    this->parameter_handler.declare_entry(
      "transport direction specification",
      "components",
      Patterns::Selection("components|2d_angle"),
      "Specification type for transport direction");
    this->parameter_handler.declare_entry(
      "azimuthal angle",
      "0.0",
      Patterns::Double(),
      "Azimuthal angle of transport direction");
    this->parameter_handler.declare_entry("polar angle",
                                          "0.0",
                                          Patterns::Double(),
                                          "Polar angle of transport direction");
    this->parameter_handler.declare_entry("transport direction x",
                                          "1.0",
                                          Patterns::Double(),
                                          "x-component of transport direction");
    this->parameter_handler.declare_entry("transport direction y",
                                          "0.0",
                                          Patterns::Double(),
                                          "y-component of transport direction");
    this->parameter_handler.declare_entry("transport direction z",
                                          "0.0",
                                          Patterns::Double(),
                                          "z-component of transport direction");
    this->parameter_handler.declare_entry(
      "normalize transport direction",
      "true",
      Patterns::Bool(),
      "Option to normalize transport direction vector");
    this->parameter_handler.declare_entry(
      "cross section", "0", Patterns::Anything(), "Cross section");
    this->parameter_handler.declare_entry(
      "source", "0", Patterns::Anything(), "Source");
    this->parameter_handler.declare_entry("source is time dependent",
                                          "false",
                                          Patterns::Bool(),
                                          "Flag that source is time-dependent");
  }
  this->parameter_handler.leave_subsection();

  // constants
  this->parameter_handler.enter_subsection("constants");
  {
    this->parameter_handler.declare_entry("incoming value",
                                          "1.0",
                                          Patterns::Double(),
                                          "incoming value for function parsers");
    this->parameter_handler.declare_entry(
      "sigma1", "0.0", Patterns::Double(), "cross section value 1");
    this->parameter_handler.declare_entry(
      "sigma2", "0.0", Patterns::Double(), "cross section value 2");
    this->parameter_handler.declare_entry(
      "x1", "0.0", Patterns::Double(), "x value 1");
    this->parameter_handler.declare_entry(
      "x2", "0.0", Patterns::Double(), "x value 2");
    this->parameter_handler.declare_entry(
      "x3", "0.0", Patterns::Double(), "x value 3");
    this->parameter_handler.declare_entry(
      "x4", "0.0", Patterns::Double(), "x value 4");
    this->parameter_handler.declare_entry(
      "y1", "0.0", Patterns::Double(), "y value 1");
    this->parameter_handler.declare_entry(
      "y2", "0.0", Patterns::Double(), "y value 2");
    this->parameter_handler.declare_entry(
      "y3", "0.0", Patterns::Double(), "y value 3");
    this->parameter_handler.declare_entry(
      "source value", "0.0", Patterns::Double(), "source for function parsers");
  }
  this->parameter_handler.leave_subsection();

  // boundary conditions
  this->parameter_handler.enter_subsection("boundary conditions");
  {
    this->parameter_handler.declare_entry(
      "dirichlet function", "0", Patterns::Anything(), "Dirichlet function");
  }
  this->parameter_handler.leave_subsection();

  // initial conditions
  this->parameter_handler.enter_subsection("initial conditions");
  {
    this->parameter_handler.declare_entry(
      "initial condition", "0", Patterns::Anything(), "Initial conditions");
  }
  this->parameter_handler.leave_subsection();

  // exact solution
  this->parameter_handler.enter_subsection("exact solution");
  {
    this->parameter_handler.declare_entry(
      "exact solution", "1", Patterns::Anything(), "Exact solution");
  }
  this->parameter_handler.leave_subsection();
}

/**
 * \brief Gets derived parameters.
 */
template <int dim>
void TransportProblemParameters<dim>::get_derived_parameters()
{
  // physics
  this->parameter_handler.enter_subsection("physics");
  {
    transport_speed = this->parameter_handler.get_double("transport speed");
    transport_direction_specification =
      this->parameter_handler.get("transport direction specification");
    azimuthal_angle = this->parameter_handler.get_double("azimuthal angle");
    polar_angle = this->parameter_handler.get_double("polar angle");
    transport_direction_x =
      this->parameter_handler.get_double("transport direction x");
    transport_direction_y =
      this->parameter_handler.get_double("transport direction y");
    transport_direction_z =
      this->parameter_handler.get_double("transport direction z");
    normalize_transport_direction =
      this->parameter_handler.get_bool("normalize transport direction");
    cross_section_string = this->parameter_handler.get("cross section");
    source_string = this->parameter_handler.get("source");
    source_is_time_dependent =
      this->parameter_handler.get_bool("source is time dependent");
  }
  this->parameter_handler.leave_subsection();

  // constants
  this->parameter_handler.enter_subsection("constants");
  {
    incoming_value = this->parameter_handler.get_double("incoming value");
    source_value = this->parameter_handler.get_double("source value");
    sigma1 = this->parameter_handler.get_double("sigma1");
    sigma2 = this->parameter_handler.get_double("sigma2");
    x1 = this->parameter_handler.get_double("x1");
    x2 = this->parameter_handler.get_double("x2");
    x3 = this->parameter_handler.get_double("x3");
    x4 = this->parameter_handler.get_double("x4");
    y1 = this->parameter_handler.get_double("y1");
    y2 = this->parameter_handler.get_double("y2");
    y3 = this->parameter_handler.get_double("y3");
  }
  this->parameter_handler.leave_subsection();

  // boundary conditions
  this->parameter_handler.enter_subsection("boundary conditions");
  {
    dirichlet_function_angularflux =
      this->parameter_handler.get("dirichlet function");
  }
  this->parameter_handler.leave_subsection();

  // initial conditions
  this->parameter_handler.enter_subsection("initial conditions");
  {
    initial_condition_angularflux =
      this->parameter_handler.get("initial condition");
  }
  this->parameter_handler.leave_subsection();

  // exact solution
  this->parameter_handler.enter_subsection("exact solution");
  {
    exact_solution_angularflux = this->parameter_handler.get("exact solution");
  }
  this->parameter_handler.leave_subsection();
}

/**
 * \brief Processes derived parameters.
 *
 * \param[in] triangulation    triangulation
 * \param[in] fe               finite element system
 * \param[in] face_quadrature  face quadrature
 */
template <int dim>
void TransportProblemParameters<dim>::process_derived_parameters(
  Triangulation<dim> & triangulation,
  const FESystem<dim> & fe,
  const QGauss<dim - 1> & face_quadrature)
{
  // constants for function parsers
  this->constants["speed"] = transport_speed;
  this->constants["incoming"] = incoming_value;
  this->constants["sigma1"] = sigma1;
  this->constants["sigma2"] = sigma2;
  this->constants["source"] = source_value;
  this->constants["x1"] = x1;
  this->constants["x2"] = x2;
  this->constants["x3"] = x3;
  this->constants["x4"] = x4;
  this->constants["y1"] = y1;
  this->constants["y2"] = y2;
  this->constants["y3"] = y3;

  // get transport direction (possibly unnormalized)
  if (transport_direction_specification == "components")
  {
    for (unsigned int d = 0; d < dim; ++d)
      if (d == 0)
        transport_direction[d] = transport_direction_x;
      else if (d == 1)
        transport_direction[d] = transport_direction_y;
      else
        transport_direction[d] = transport_direction_z;
  }
  else if (transport_direction_specification == "2d_angle")
  {
    // assert not 3-D
    Assert(dim != 3, ExcImpossibleInDim(dim));

    // compute direction components
    const double azimuthal_angle_radians = azimuthal_angle * numbers::PI / 180.0;
    transport_direction[0] = std::cos(azimuthal_angle_radians);
    if (dim == 2)
      transport_direction[1] = std::sin(azimuthal_angle_radians);
  }
  else
  {
    Assert(false, ExcNotImplemented());
  }

  // normalize transport direction
  double transport_direction_magnitude = 0.0;
  for (unsigned int d = 0; d < dim; ++d)
    transport_direction_magnitude += std::pow(transport_direction[d], 2);
  transport_direction_magnitude = std::sqrt(transport_direction_magnitude);
  for (unsigned int d = 0; d < dim; ++d)
    transport_direction[d] /= transport_direction_magnitude;

  // initialize cross section function
  cross_section_function.initialize(FunctionParser<dim>::default_variable_names(),
                                    cross_section_string,
                                    this->constants,
                                    false);

  // initialize source function
  source_function.initialize(FunctionParser<dim>::default_variable_names() + ",t",
                             source_string,
                             this->constants,
                             true);

  // store Dirichlet BC strings
  this->dirichlet_function_strings.resize(1);
  this->dirichlet_function_strings[0] = dirichlet_function_angularflux;

  // store initial condition strings
  this->initial_conditions_strings.resize(1);
  this->initial_conditions_strings[0] = initial_condition_angularflux;

  // store exact solution strings
  this->exact_solution_strings.resize(1);
  this->exact_solution_strings[0] = exact_solution_angularflux;

  // set boundary IDs if using the incoming boundary specification
  if (this->boundary_id_scheme == "incoming")
    set_boundary_ids_incoming(triangulation, fe, face_quadrature);
}

/**
 * \brief Sets the boundary IDs for incoming boundary specification.
 *
 * \pre \c FEFaceValues object must use the flag \c update_normal_vectors.
 *
 * \param[in] triangulation    triangulation
 * \param[in] fe               finite element system
 * \param[in] face_quadrature  face quadrature
 */
template <int dim>
void TransportProblemParameters<dim>::set_boundary_ids_incoming(
  Triangulation<dim> & triangulation,
  const FESystem<dim> & fe,
  const QGauss<dim - 1> & face_quadrature)
{
  // compute faces per cell
  const unsigned int faces_per_cell = GeometryInfo<dim>::faces_per_cell;

  // reset boundary indicators to one
  typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
  typename Triangulation<dim>::cell_iterator endc = triangulation.end();
  for (; cell != endc; ++cell)
    for (unsigned int face = 0; face < faces_per_cell; ++face)
      if (cell->face(face)->at_boundary())
        cell->face(face)->set_boundary_id(1);

  // create FE face values for evaluating normal vectors
  FEFaceValues<dim> fe_face_values(fe, face_quadrature, update_normal_vectors);

  // loop over cells
  cell = triangulation.begin();
  endc = triangulation.end();
  for (; cell != endc; ++cell)
  {
    // loop over faces of cell
    for (unsigned int face = 0; face < faces_per_cell; ++face)
    {
      // if face is at boundary
      if (cell->face(face)->at_boundary())
      {
        // reinitialize FE face values
        fe_face_values.reinit(cell, face);

        // determine if the transport flux is incoming through this face;
        //  it isn't necessary to loop over all face quadrature points because
        //  the transport direction and normal vector are the same at each
        //  quadrature point; therefore, quadrature point 0 is arbitrarily
        //  chosen
        const double small = -1.0e-12;
        const double n_dot_omega = fe_face_values.normal_vector(0) * transport_direction;
        if (n_dot_omega < small)
        {
          // mark boundary as incoming flux boundary: indicator 0
          cell->face(face)->set_boundary_id(0);
        }
      }
    }
  }
}
