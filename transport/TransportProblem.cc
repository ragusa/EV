/** \brief constructor for TransportProblem class
 */
template<int dim>
TransportProblem<dim>::TransportProblem(const TransportParameters &parameters) :
      parameters(parameters),
      dof_handler(triangulation),
      degree(parameters.degree),
      fe(FE_Q<dim>(degree), 1),
      flux(0),
      dofs_per_cell(fe.dofs_per_cell),
      faces_per_cell(GeometryInfo<dim>::faces_per_cell),
      cell_quadrature(parameters.n_quadrature_points),
      face_quadrature(parameters.n_quadrature_points),
      n_q_points_cell(cell_quadrature.size()),
      n_q_points_face(face_quadrature.size()),
      transport_direction(0.0),
      incoming_boundary(1),
      has_exact_solution(false),
      incoming_flux_value(0.0),
      domain_volume(0.0)
{
}

/** \brief destructor for TransportProblem class
 */
template<int dim>
TransportProblem<dim>::~TransportProblem() {
   dof_handler.clear();
}

/** \brief initialize system
 */
template<int dim>
void TransportProblem<dim>::initialize_system()
{
   // process problem ID
   process_problem_ID();

   // decide that transport direction is the unit x vector
   transport_direction[0] = 1.0;

   // determine function variables based on dimension;
   // this is to be used in initialization of function parser objects
   std::string variables;
   if (dim == 1)
      variables = "x";
   else if (dim == 2)
      variables = "x,y";
   else if (dim == 3)
      variables = "x,y,z";
   else
      Assert(false,ExcInvalidState());

   // add time variable if not steady state
   if (!parameters.is_steady_state)
      variables += ",t";

   // create function parser constants
   function_parser_constants["incoming"] = incoming_flux_value;
   function_parser_constants["pi"]       = numbers::PI;
   function_parser_constants["x_min"]    = x_min;
   function_parser_constants["x_mid"]    = x_min + 0.5*(x_max-x_min);

   // initialize exact solution function
   if (has_exact_solution)
      exact_solution.initialize(variables,
                                exact_solution_string,
                                function_parser_constants,
                                true);

   // initialize initial conditions function
   if (!parameters.is_steady_state)
      initial_conditions.initialize(variables,
                                    initial_conditions_string,
                                    function_parser_constants,
                                    true);

   // initialize source function
   source_function.initialize(variables,
                              source_string,
                              function_parser_constants,
                              true);

   // initialize cross section function
   cross_section_function.initialize(variables,
                                     cross_section_string,
                                     function_parser_constants,
                                     true);

   // create grid for initial refinement level
   GridGenerator::hyper_cube(triangulation, x_min, x_max);
   domain_volume = std::pow((x_max-x_min),dim);
   triangulation.refine_global(parameters.initial_refinement_level);
   n_cells = triangulation.n_active_cells();

   // get initial dt size - CFL imposed in setup_system()
   dt_nominal = parameters.time_step_size;
}

/** \brief process problem ID
 */
template<int dim>
void TransportProblem<dim>::process_problem_ID()
{
   switch (parameters.problem_id)
   {
      case 1: { // incoming flux, constant cross section, constant source
         Assert(dim < 3,ExcNotImplemented());

         x_min = 0.0;
         x_max = 10.0;

         incoming_flux_value = 1.0;
         function_parser_constants["incoming"]  = 1.0;

         cross_section_string = "1.0";
         function_parser_constants["sigma"]  = 1.0;

         source_time_dependent = false;
         source_string = "0";
         function_parser_constants["source"] = 0.0;

         if (parameters.is_steady_state) {
            has_exact_solution = true;
            exact_solution_string = "source/sigma + (incoming - source/sigma)*exp(-sigma*(x-x_min))";
         } else {
            has_exact_solution = false;
         }

         initial_conditions_string = "0";

         break;
      } case 2: { // high cross section in half (1-D) or quarter (2-D) of domain, no source, incoming flux
         Assert(dim < 3,ExcNotImplemented());

         x_min = 0.0;
         x_max = 10.0;

         incoming_flux_value = 1.0;
         function_parser_constants["incoming"]  = 1.0;

         cross_section_string = "if(x<x_mid, 0, 1)";
         function_parser_constants["sigma"]  = 1.0;

         source_time_dependent = false;
         source_string = "0";
         function_parser_constants["source"] = 0.0;

         has_exact_solution = true;
         if (parameters.is_steady_state) { // steady-state
            if (dim == 1) // 1-D
               exact_solution_string = "if(x<x_mid, incoming, exp(-sigma*(x-x_mid)))";
            else if (dim == 2) // 2-D
               exact_solution_string = "if(y<x_mid, incoming, if(x<x_mid, incoming, exp(-sigma*(x-x_mid))))";
            else
               Assert(false,ExcNotImplemented());
         } else { // transient
            if (dim == 1) // 1-D
               exact_solution_string = "if(x-t<x_min, if(x<x_mid, incoming, exp(-sigma*(x-x_mid))), 0)";
            else if (dim == 2) // 2-D
               exact_solution_string = "if(x-t<x_min, if(y<x_mid, incoming, if(x<x_mid, incoming, exp(-sigma*(x-x_mid)))), 0)";
            else
               Assert(false,ExcNotImplemented());
         }
         initial_conditions_string = "0";
         break;
      } case 3: { // 1-D linear advection problem with step IC
         Assert(dim == 1,ExcNotImplemented());

         x_min = 0.0;
         x_max = 1.0;
         cross_section_string = "0";
         source_string = "0";
         source_time_dependent = false;
         incoming_flux_value = 1.0;
         has_exact_solution = true;
         if (parameters.is_steady_state)
            exact_solution_string = "1.0";
         else
            exact_solution_string = "if((x-t)<0,incoming,0)";
         initial_conditions_string = "if(x<0,incoming,0)";
         break;
      } case 5: { // MMS-1
         Assert(dim == 1,ExcNotImplemented()); // assume 1-D

         x_min = 0.0;
         x_max = 1.0;
         incoming_flux_value = 0.0;
         cross_section_string = "1.0";
         has_exact_solution = true;
         // use different MS for steady-state vs. transient problem
         if (parameters.is_steady_state) {
            exact_solution_string = "sin(pi*x)"; // assume omega_x = 1 and c = 1
            source_string = "pi*cos(pi*x) + sin(pi*x)";
            source_time_dependent = false;
         } else {
            exact_solution_string = "t*sin(pi*x)"; // assume omega_x = 1 and c = 1
            source_string = "sin(pi*x) + pi*t*cos(pi*x) + t*sin(pi*x)";
            source_time_dependent = true;
            initial_conditions_string = exact_solution_string;
         }
         break;
      } case 6: { // MMS-2
         Assert(dim == 1,ExcNotImplemented()); // assume 1-D
         Assert(not parameters.is_steady_state,ExcNotImplemented()); // assume not steady-state

         x_min = 0.0;
         x_max = 1.0;
         incoming_flux_value = 0.0;
         cross_section_string = "1";
         has_exact_solution = true;
         exact_solution_string = "x*exp(-t)"; // assume omega_x = 1 and c = 1
         source_string = "-x*exp(-t) + exp(-t) + x*exp(-t)";
         source_time_dependent = true;
         initial_conditions_string = exact_solution_string;
         break;
      } case 7: { // MMS-3
         Assert(dim == 1,ExcNotImplemented()); // assume 1-D
         Assert(not parameters.is_steady_state,ExcNotImplemented()); // assume not steady-state

         x_min = 0.0;
         x_max = 1.0;
         incoming_flux_value = 0.0;
         cross_section_string = "1";
         has_exact_solution = true;
         exact_solution_string = "t^4*sin(pi*x)"; // assume omega_x = 1 and c = 1
         source_string = "4*t^3*sin(pi*x) + pi*t^4*cos(pi*x) + t^4*sin(pi*x)";
         source_time_dependent = true;
         initial_conditions_string = exact_solution_string;
         break;
      } default: {
         Assert(false,ExcNotImplemented());
         break;
      }
   }
}

/** \brief set up the problem before assembly of the linear system
 */
template<int dim>
void TransportProblem<dim>::setup_system()
{
   // distribute dofs
   dof_handler.distribute_dofs(fe);
   // get number of dofs
   n_dofs = dof_handler.n_dofs();

   // set boundary indicators to distinguish incoming boundary
   set_boundary_indicators();
   // determine Dirichlet nodes
   get_dirichlet_nodes();

   // compute minimum cell diameter for CFL condition
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   // reset minimum cell diameter to an arbitrary cell diameter such as the first one
   minimum_cell_diameter = cell->diameter();
   // loop over cells to determine minimum cell diameter
   for (; cell != endc; ++cell)
      minimum_cell_diameter = std::min(minimum_cell_diameter, cell->diameter());

   // clear constraint matrix and make hanging node constraints
   constraints.clear();
   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
   constraints.close();

   // create sparsity pattern for system matrix and mass matrices
   CompressedSparsityPattern compressed_constrained_sparsity_pattern(n_dofs);
   DoFTools::make_sparsity_pattern(dof_handler,
                                   compressed_constrained_sparsity_pattern,
                                   constraints,
                                   false);
   constrained_sparsity_pattern.copy_from(compressed_constrained_sparsity_pattern);

   // reinitialize system matrix and mass matrices
   system_matrix            .reinit(constrained_sparsity_pattern);
   low_order_ss_matrix      .reinit(constrained_sparsity_pattern);
   high_order_ss_matrix     .reinit(constrained_sparsity_pattern);
   inviscid_ss_matrix       .reinit(constrained_sparsity_pattern);
   low_order_viscous_matrix .reinit(constrained_sparsity_pattern);
   consistent_mass_matrix   .reinit(constrained_sparsity_pattern);
   lumped_mass_matrix       .reinit(constrained_sparsity_pattern);
   if ((parameters.scheme_option == 2)||(parameters.scheme_option == 3))
   {
      high_order_viscous_matrix.reinit(constrained_sparsity_pattern);
      flux_correction_matrix.reinit(constrained_sparsity_pattern);
   }
      
   // reinitialize auxiliary matrices
   CompressedSparsityPattern compressed_unconstrained_sparsity_pattern(n_dofs);
   DoFTools::make_sparsity_pattern(dof_handler, compressed_unconstrained_sparsity_pattern);
   unconstrained_sparsity_pattern.copy_from(compressed_unconstrained_sparsity_pattern);
   viscous_bilinear_forms            .reinit(unconstrained_sparsity_pattern);
   // compute viscous bilinear forms
   compute_viscous_bilinear_forms();

   // reinitialize solution vector, system matrix, and rhs
   old_solution.reinit(n_dofs);
   new_solution.reinit(n_dofs);
   system_rhs  .reinit(n_dofs);
   ss_rhs      .reinit(n_dofs);
   R_plus      .reinit(n_dofs);
   R_minus     .reinit(n_dofs);
   Q_plus      .reinit(n_dofs);
   Q_minus     .reinit(n_dofs);
   flux_correction_vector.reinit(n_dofs);

   // reinitialize max principle bounds
   min_values.reinit(n_dofs);
   max_values.reinit(n_dofs);

   // reinitialize viscosities
   entropy_viscosity   .reinit(n_cells);
   low_order_viscosity .reinit(n_cells);
   high_order_viscosity.reinit(n_cells);

   // assemble mass matrices
   assemble_mass_matrices();

   // compute inviscid system matrix and steady-state right hand side (ss_rhs)
   assemble_inviscid_ss_matrix();
   if (not source_time_dependent)
      assemble_ss_rhs(0.0);
      
   // compute low-order viscous matrix and low-order steady-state matrix
   compute_low_order_viscosity();
   compute_viscous_matrix(low_order_viscosity,low_order_viscous_matrix);
   low_order_ss_matrix.copy_from(inviscid_ss_matrix);
   low_order_ss_matrix.add(1.0,low_order_viscous_matrix);

   // enforce CFL condition on nominal dt size
   CFL_nominal = enforce_CFL_condition(dt_nominal);
}

/** \brief Assembles the mass matrix, either consistent or lumped.
 */
template <int dim>
void TransportProblem<dim>::assemble_mass_matrices()
{
   // assemble lumped mass matrix
   //----------------------------
   FEValues<dim> fe_values (fe, cell_quadrature, update_values | update_JxW_values);

   std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
   FullMatrix<double> local_mass (dofs_per_cell, dofs_per_cell);

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell)
   {
      fe_values.reinit(cell);
      cell->get_dof_indices(local_dof_indices);

      local_mass = 0.0;

      // compute local contribution
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
               local_mass(i,i) +=  fe_values.shape_value(i,q)
                                  *fe_values.shape_value(j,q)
                                  *fe_values.JxW(q);
            }

      // add to global mass matrix with contraints
      constraints.distribute_local_to_global (local_mass, local_dof_indices, lumped_mass_matrix);
   }

   // assemble consistent mass matrix
   //----------------------------
   Function<dim> *dummy_function = 0;
   MatrixTools::create_mass_matrix(dof_handler,
                                   cell_quadrature,
                                   consistent_mass_matrix,
                                   dummy_function,
                                   constraints);
}

/** \brief Computes viscous bilinear forms, to be used in the computation of
 *         maximum-principle preserving low order viscosity.
 *
 *         Each element of the resulting matrix, \f$B_{i,j}\f$ is computed as
 *         follows:
 *         \f[
 *            B_{i,j} = \sum_{K\subset S_{i,j}}b_K(\varphi_i,\varphi_j)
 *         \f]
 */
template <int dim>
void TransportProblem<dim>::compute_viscous_bilinear_forms()
{
   viscous_bilinear_forms = 0; // zero out matrix

   // dofs per cell
   unsigned int dofs_per_cell = fe.dofs_per_cell;

   // local dof indices
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // query cell volume
      double cell_volume = cell->measure();

      // add local bilinear forms to global matrix
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
         for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            double b_cell = 0.0;
            if (j == i) {
               b_cell = cell_volume;
            } else {
               b_cell = -1.0/(dofs_per_cell-1.0)*cell_volume;
            }
            viscous_bilinear_forms.add(local_dof_indices[i],
                                       local_dof_indices[j],
                                       b_cell);
         }
      }
   }
}

/** \brief Assemble the inviscid system matrix. The inviscid steady-state matrix
 *         is independent of the solution and independent of time and thus needs
 *         to be called only once per level of mesh refinement.
 */
template<int dim>
void TransportProblem<dim>::assemble_inviscid_ss_matrix()
{
   inviscid_ss_matrix = 0;

   // FE values, for assembly terms
   FEValues<dim> fe_values(fe, cell_quadrature,
         update_values | update_gradients | update_quadrature_points
               | update_JxW_values);

   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // total cross section values at each quadrature point on cell
   std::vector<double> total_cross_section_values(n_q_points_cell);

   // cell iterator
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   // loop over cells
   unsigned int i_cell = 0;
   for (cell = dof_handler.begin_active(); cell != endc; ++cell, ++i_cell) {
      // initialize local matrix and rhs to zero
      cell_matrix = 0;

      // reinitialize FE values
      fe_values.reinit(cell);

      // get quadrature points on cell
      std::vector<Point<dim> > points(n_q_points_cell);
      points = fe_values.get_quadrature_points();

      // get cross section values for all quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         total_cross_section_values[q] = cross_section_function.value(points[q]);

      // compute cell contributions to global system
      // ------------------------------------------------------------------
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
               // store integrals of divergence and total interaction term
               // so that they may be used in computation of max-principle
               // preserving viscosity
               cell_matrix(i,j) += (
                  // divergence term
                  fe_values[flux].value(i, q)
                     * transport_direction
                     * fe_values[flux].gradient(j, q) +
                  // total interaction term
                  fe_values[flux].value(i, q)
                     * total_cross_section_values[q]
                     * fe_values[flux].value(j, q)
                  ) * fe_values.JxW(q);
            } // end j
         } // end i
      } // end q

      // aggregate local matrix and rhs to global matrix and rhs
      constraints.distribute_local_to_global(cell_matrix,
                                             local_dof_indices,
                                             inviscid_ss_matrix);
   } // end cell
} // end assembly

/** \brief Assemble the steady-state rhs.
 */
template<int dim>
void TransportProblem<dim>::assemble_ss_rhs(const double &t)
{
   // reset steady-state rhs
   ss_rhs = 0;

   // set the time to be used in the source function
   source_function.set_time(t);

   // FE values, for assembly terms
   FEValues<dim> fe_values(fe, cell_quadrature,
         update_values | update_gradients | update_quadrature_points
               | update_JxW_values);

   Vector<double>            cell_rhs   (dofs_per_cell);

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // source values at each quadrature point on cell
   std::vector<Point<dim> > points(n_q_points_cell);
   std::vector<double> source_values(n_q_points_cell);

   // cell iterator
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   // loop over cells
   for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
      // initialize local rhs to zero
      cell_rhs = 0;

      // reinitialize FE values
      fe_values.reinit(cell);

      // get quadrature points on cell
      points = fe_values.get_quadrature_points();

      // get total source for all quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         source_values[q] = source_function.value(points[q]);

      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
            cell_rhs(i) +=
               fe_values[flux].value(i, q)
                  * source_values[q] * fe_values.JxW(q);

      // aggregate local matrix and rhs to global matrix and rhs
      constraints.distribute_local_to_global(cell_rhs, local_dof_indices, ss_rhs);
   } // end cell
} // end assembly

/** \brief Applies Dirichlet boundary conditions to a linear system A*x=b.
 *  \param [in,out] A linear system matrix
 *  \param [in,out] x linear system solution vector
 *  \param [in,out] b linear system rhs vector
 */
template <int dim>
void TransportProblem<dim>::apply_Dirichlet_BC(SparseMatrix<double> &A,
                                               Vector<double>       &x,
                                               Vector<double>       &b)
{
   // apply Dirichlet boundary condition
   std::map<unsigned int, double> boundary_values;
   VectorTools::interpolate_boundary_values(dof_handler,
                                            incoming_boundary,
                                            ConstantFunction<dim>(incoming_flux_value, 1),
                                            boundary_values);
   MatrixTools::apply_boundary_values(boundary_values,
                                      A,
                                      x,
                                      b);
}

/** \brief adds the viscous bilinear form for the maximum-principle preserving viscosity
 */
template <int dim>
void TransportProblem<dim>::compute_viscous_matrix(const Vector<double> &viscosity,
                                                   SparseMatrix<double> &viscous_matrix)
{
   // reset viscous matrix
   viscous_matrix = 0;

   unsigned int i_cell = 0;
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell, ++i_cell) {
      // reset cell matrix to zero
      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      cell_matrix = 0;

      // compute cell volume
      double cell_volume = cell->measure();
      // compute cell contribution to global system matrix
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
         for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            double viscous_bilinear_form;
            if (j == i)
               viscous_bilinear_form = cell_volume;
            else
               viscous_bilinear_form = cell_volume/(1.0 - dofs_per_cell);
 
            cell_matrix(i,j) += viscosity(i_cell) * viscous_bilinear_form;
         }
      }

      // aggregate local matrix and rhs to global matrix and rhs
      std::vector<unsigned int> local_dof_indices(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix,
                                             local_dof_indices,
                                             viscous_matrix);

   }
}

/** \brief Sets the boundary indicators for each boundary face.
 *
 *         The Dirichlet BC is applied only to the incoming boundary, so the transport
 *         direction is compared against the normal vector of the face.
 */
template <int dim>
void TransportProblem<dim>::set_boundary_indicators()
{
   // cell iterator
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   // reset boundary indicators to zero
   for (cell = dof_handler.begin_active(); cell != endc; ++cell)
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
         if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_indicator(0);

   // FE face values
   FEFaceValues<dim> fe_face_values(fe, face_quadrature, update_normal_vectors);
   // loop over cells
   for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
      // loop over faces of cell
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face) {
         // if face is at boundary
         if (cell->face(face)->at_boundary()) {
            // reinitialize FE face values
            fe_face_values.reinit(cell, face);
            // determine if the transport flux is incoming through this face;
            //  it isn't necessary to loop over all face quadrature points because
            //  the transport direction and normal vector are the same at each
            //  quadrature point; therefore, quadrature point 0 is arbitrarily chosen
            double small = -1.0e-12;
            if (fe_face_values.normal_vector(0) * transport_direction < small) {
               // mark boundary as incoming flux boundary: indicator 1
               cell->face(face)->set_boundary_indicator(incoming_boundary);
            }
         }
      }
   }
}

/** \brief Computes the low-order viscosity for each cell.
 */
template <int dim>
void TransportProblem<dim>::compute_low_order_viscosity()
{
   // local dof indices
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   // loop over cells to compute low-order viscosity at each quadrature point
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   unsigned int i_cell = 0;
   for (; cell != endc; ++cell, ++i_cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // initialize to zero for max() function
      double low_order_viscosity_cell = 0.0;
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
         for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            if (i != j) {
               low_order_viscosity_cell = std::max(low_order_viscosity_cell,
                  std::max(0.0,inviscid_ss_matrix(local_dof_indices[i],local_dof_indices[j]))/
                  (-viscous_bilinear_forms(local_dof_indices[i],local_dof_indices[j])));
            }
         } 
      }
      low_order_viscosity(i_cell) = low_order_viscosity_cell;
   }
}

/** \brief refine the grid
 */
template<int dim>
void TransportProblem<dim>::refine_grid() {
   // refine adaptively or globally
   if (parameters.use_adaptive_mesh_refinement) {
      Vector<float> estimated_error_per_cell(n_cells);

      KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim - 1>(3),
            typename FunctionMap<dim>::type(), new_solution,
            estimated_error_per_cell);

      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
            estimated_error_per_cell, 0.3, 0.03);

      triangulation.execute_coarsening_and_refinement();
   } else
      triangulation.refine_global(1);

   // update number of cells
   n_cells = triangulation.n_active_cells();
}

/** \brief output the grid of the given cycle
 */
template<int dim>
void TransportProblem<dim>::output_grid() const
{
   // create output filestream
   std::string filename = "grid.eps";
   std::ofstream output(filename.c_str());

   // write grid to eps
   GridOut grid_out;
   grid_out.write_eps(triangulation, output);
}

/** \brief run the problem
 */
template<int dim>
void TransportProblem<dim>::run()
{
   // initialize problem
   initialize_system();

   // loop over refinement cycles
   for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; ++cycle)
   {
      // refine
      if (cycle != 0) {
         // refine mesh if user selected the option
         if (parameters.refinement_option != 2) // "2" corresponds to "time refinement only"
            refine_grid();

         // refine time if user selected the option
         if (parameters.refinement_option != 1) // "1" corresponds to "space refinement only"
            dt_nominal = dt_nominal * parameters.time_refinement_factor;
      }

      // setup system - distribute finite elements, reintialize matrices and vectors
      setup_system();

      // print information
      std::cout << std::endl << "Cycle " << cycle << ':' << std::endl;
      std::cout << "   Number of active cells:       " << n_cells << std::endl;
      std::cout << "   Number of degrees of freedom: " << n_dofs  << std::endl;
      if (not parameters.is_steady_state) {
         std::cout << "   Nominal time step size: " << dt_nominal << std::endl;
         std::cout << "   Nominal CFL number: " << CFL_nominal << std::endl;
      }

      // create entropy viscosity object
      EntropyViscosity<dim> EV(fe,
                               n_cells,
                               dof_handler,
                               cell_quadrature,
                               face_quadrature,
                               transport_direction,
                               cross_section_function,
                               source_function,
                               parameters.entropy_string,
                               parameters.entropy_derivative_string,
                               parameters.entropy_residual_coefficient,
                               parameters.jump_coefficient,
                               domain_volume);                      

      // create linear solver object
      LinearSolver linear_solver(parameters.linear_solver_option, constraints);

      // if problem is steady-state, then just do one solve; else loop over time
      if (parameters.is_steady_state) { // run steady-state problem
         solve_steady_state();
      } else {                          // run transient problem
         // interpolate initial conditions
         initial_conditions.set_time(0.0);
         VectorTools::interpolate(dof_handler,
                                  initial_conditions,
                                  new_solution);
         constraints.distribute(new_solution);

         // if last cycle, output initial conditions if user requested
         if ((cycle == parameters.n_refinement_cycles-1) and
            (parameters.output_initial_solution))
            output_solution(new_solution,
                            dof_handler,
                            "initial_solution",
                            false);

         // set old solution to the current solution
         old_solution = new_solution;

         // time loop
         double t_new = 0.0;
         double t_old = 0.0;
         const double t_end = parameters.end_time;
         bool in_transient = true;
         bool DMP_satisfied_at_all_steps = true;
         double old_dt = dt_nominal; // time step size of previous time step, needed for entropy viscosity
         Vector<double> tmp_vector(n_dofs);
         unsigned int n = 1; // time step index
         while (in_transient)
         {
            // compute high-order viscosity and viscous matrix if needed
            if ((parameters.scheme_option == 2)||(parameters.scheme_option == 3)) {
               if (n == 1) {
                  compute_viscous_matrix(low_order_viscosity,high_order_viscous_matrix);
               } else {
                  entropy_viscosity = EV.compute_entropy_viscosity(new_solution,old_solution,old_dt);
                  for (unsigned int i_cell = 0; i_cell < n_cells; ++i_cell)
                     high_order_viscosity(i_cell) = std::min(entropy_viscosity(i_cell),low_order_viscosity(i_cell));
                  compute_viscous_matrix(high_order_viscosity,high_order_viscous_matrix);
               }
            }

            // update old solution to previous step's new solution
            // Note that this is done here because the old old_solution is needed in the
            // computation of the entropy viscosity time derivative term, which occurs in
            // compute_entropy_viscosity(), so the following update to old_solution needs to come
            // sometime after compute_entropy_viscosity() but before old_solution is used in the
            // assembly of the transient residual.
            old_solution = new_solution;

            // shorten dt if new time would overshoot end time, then update t_new
            double dt = dt_nominal;
            if (t_old + dt >= t_end) {
               dt = t_end - t_old;
               in_transient = false;
            }
            t_new = t_old + dt;
            std::cout << "   time step " << n << ": t = " << t_old << "->" << t_new << std::endl;

            switch (parameters.scheme_option)
            {
               case 0: { // unmodified Galerkin scheme
                  switch (parameters.time_integrator_option)
                  {
                     case 1: { // forward Euler
                        transient_step(old_solution,old_solution,consistent_mass_matrix,inviscid_ss_matrix,dt,t_old,linear_solver);

                        break;
                     } case 2: { // SSPRK22
                        // stage 1
                        tmp_vector = old_solution;
                        transient_step(tmp_vector,tmp_vector,consistent_mass_matrix,inviscid_ss_matrix,dt,t_old,linear_solver);
                        // stage 2
                        tmp_vector = new_solution;
                        transient_step(tmp_vector,tmp_vector,consistent_mass_matrix,inviscid_ss_matrix,dt,t_new,linear_solver);
                        // final combination
                        tmp_vector = new_solution;
                        new_solution = 0;
                        new_solution.add(0.5,old_solution,0.5,tmp_vector);
                        constraints.distribute(new_solution);

                        break;
                     } case 3: { // SSPRK33
                        // stage 1
                        tmp_vector = old_solution;
                        transient_step(tmp_vector,tmp_vector,consistent_mass_matrix,inviscid_ss_matrix,dt,t_old,linear_solver);
                        // stage 2
                        tmp_vector = new_solution;
                        transient_step(tmp_vector,tmp_vector,consistent_mass_matrix,inviscid_ss_matrix,dt,t_new,linear_solver);
                        // stage 3
                        tmp_vector = 0;
                        tmp_vector.add(0.75,old_solution,0.25,new_solution);
                        transient_step(tmp_vector,tmp_vector,consistent_mass_matrix,inviscid_ss_matrix,dt,t_old+0.5*dt,linear_solver);
                        // final combination
                        tmp_vector = new_solution;
                        new_solution = 0;
                        new_solution.add(1./3.,old_solution,2./3.,tmp_vector);
                        constraints.distribute(new_solution);

                        break;
                     } default: {
                        Assert(false, ExcNotImplemented());
                        break;
                     }
                  }

                  break;
               } case 1: { // solve low-order system
                  switch (parameters.time_integrator_option)
                  {
                     case 1: { // forward Euler
                        transient_step(old_solution,old_solution,lumped_mass_matrix,low_order_ss_matrix,dt,t_old,linear_solver);

                        break;
                     } case 2: { // SSPRK22
                        // stage 1
                        tmp_vector = old_solution;
                        transient_step(tmp_vector,tmp_vector,lumped_mass_matrix,low_order_ss_matrix,dt,t_old,linear_solver);
                        // stage 2
                        tmp_vector = new_solution;
                        transient_step(tmp_vector,tmp_vector,lumped_mass_matrix,low_order_ss_matrix,dt,t_new,linear_solver);
                        // final combination
                        tmp_vector = new_solution;
                        new_solution = 0;
                        new_solution.add(0.5,old_solution,0.5,tmp_vector);
                        constraints.distribute(new_solution);

                        break;
                     } case 3: { // SSPRK33
                        // stage 1
                        tmp_vector = old_solution;
                        transient_step(tmp_vector,tmp_vector,lumped_mass_matrix,low_order_ss_matrix,dt,t_old,linear_solver);
                        // stage 2
                        tmp_vector = new_solution;
                        transient_step(tmp_vector,tmp_vector,lumped_mass_matrix,low_order_ss_matrix,dt,t_new,linear_solver);
                        // stage 3
                        tmp_vector = 0;
                        tmp_vector.add(0.75,old_solution,0.25,new_solution);
                        transient_step(tmp_vector,tmp_vector,lumped_mass_matrix,low_order_ss_matrix,dt,t_old+0.5*dt,linear_solver);
                        // final combination
                        tmp_vector = new_solution;
                        new_solution = 0;
                        new_solution.add(1./3.,old_solution,2./3.,tmp_vector);
                        constraints.distribute(new_solution);

                        break;
                     } default: {
                        Assert(false, ExcNotImplemented());
                        break;
                     }
                  }
                  
/*
                  // compute max principle min and max values
                  compute_max_principle_bounds(dt);
                  // check that local discrete maximum principle is satisfied at all time steps
                  bool DMP_satisfied_this_step = check_max_principle(dt,false);
                  DMP_satisfied_at_all_steps = DMP_satisfied_at_all_steps and DMP_satisfied_this_step;
*/

                  break;
               } case 2: { // high-order system with entropy viscosity
                  switch (parameters.time_integrator_option)
                  {
                     case 1: { // forward Euler
                        high_order_ss_matrix.copy_from(inviscid_ss_matrix);
                        high_order_ss_matrix.add(1.0,high_order_viscous_matrix);
                        transient_step(old_solution,old_solution,consistent_mass_matrix,high_order_ss_matrix,dt,t_old,linear_solver);

                        break;
                     } case 2: { // SSPRK33
                        Assert(false, ExcNotImplemented());
                        break;
                     } default: {
                        Assert(false, ExcNotImplemented());
                        break;
                     }
                  }

                  break;
               } case 3: { // FCT
                  Assert(false, ExcNotImplemented());
/*
                  // SSPRK33 not yet implemented
                  if (parameters.time_integrator_option != 1)
                     Assert(false, ExcNotImplemented());

                  // solve high-order system
                  // ------------------------------------------------------------
                  // form transient rhs: system_rhs = M*u_old + dt*(ss_rhs - A*u_old)
                  system_rhs = 0;
                  system_rhs.add(dt, ss_rhs); //       now, system_rhs = dt*(ss_rhs)
                  consistent_mass_matrix.vmult(tmp_vector, old_solution);
                  system_rhs.add(1.0, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs)
                  inviscid_ss_matrix.vmult(tmp_vector, old_solution);
                  system_rhs.add(-dt, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs - A*u_old)
      
                  // add viscous matrix to rhs
                  high_order_viscous_matrix.vmult(tmp_vector, old_solution);
                  system_rhs.add(-dt, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs - (A+D)*u_old)

                  // solve the linear system M*u_new = system_rhs
                  system_matrix.copy_from(consistent_mass_matrix);
                  apply_Dirichlet_BC(system_matrix, new_solution, system_rhs); // apply Dirichlet BC
                  solve_linear_system(system_matrix, system_rhs);
            
                  // distribute constraints
                  constraints.distribute(new_solution);

                  // solve FCT system
                  // ------------------------------------------------------------
                  // form transient rhs: system_rhs = M*u_old + dt*(ss_rhs - (A+D)*u_old)
                  system_rhs = 0;
                  system_rhs.add(dt, ss_rhs); //       now, system_rhs = dt*(ss_rhs)
                  lumped_mass_matrix.vmult(tmp_vector, old_solution);
                  system_rhs.add(1.0, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs)
                  low_order_ss_matrix.vmult(tmp_vector, old_solution);
                  system_rhs.add(-dt, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs - A*u_old)

                  // add limited flux sum
                  assemble_flux_correction_matrix(dt);

                  // compute max principle min and max values
                  compute_max_principle_bounds(dt);

                  // compute Q+
                  Q_plus = 0;
                  lumped_mass_matrix.vmult(tmp_vector, max_values);
                  Q_plus.add(1.0/dt, tmp_vector);
                  lumped_mass_matrix.vmult(tmp_vector, old_solution);
                  Q_plus.add(-1.0/dt,tmp_vector);
                  low_order_ss_matrix.vmult(tmp_vector, old_solution);
                  Q_plus.add(1.0,tmp_vector);
                  Q_plus.add(-1.0,ss_rhs);
                  
                  // compute Q-
                  Q_minus = 0;
                  lumped_mass_matrix.vmult(tmp_vector, max_values);
                  Q_minus.add(1.0/dt, tmp_vector);
                  lumped_mass_matrix.vmult(tmp_vector, old_solution);
                  Q_minus.add(-1.0/dt,tmp_vector);
                  low_order_ss_matrix.vmult(tmp_vector, old_solution);
                  Q_minus.add(1.0,tmp_vector);
                  Q_minus.add(-1.0,ss_rhs);

                  // compute limited flux correction sum and add it to rhs
                  compute_limiting_coefficients();
                  system_rhs.add(dt, flux_correction_vector);   // now, system_rhs = M*u_old + dt*(ss_rhs - A*u_old + f)

                  // solve the linear system M*u_new = system_rhs
                  system_matrix.copy_from(lumped_mass_matrix);
                  apply_Dirichlet_BC(system_matrix, new_solution, system_rhs); // apply Dirichlet BC
                  solve_linear_system(system_matrix, system_rhs);

                  // distribute constraints
                  constraints.distribute(new_solution);

                  // check that local discrete maximum principle is satisfied at all time steps
                  bool DMP_satisfied_this_step = check_max_principle(dt,false);
                  DMP_satisfied_at_all_steps = DMP_satisfied_at_all_steps and DMP_satisfied_this_step;
*/

                  break;
               } default: {
                  Assert(false, ExcNotImplemented());
                  break;
               }

            }

            // update old time and time step size
            t_old = t_new;
            old_dt = dt;

            // increment time step index
            n++;
         }

         // report if DMP was satisfied at all time steps
         if (parameters.scheme_option == 3)
         {
            // report if local discrete maximum principle is satisfied at all time steps
            if (DMP_satisfied_at_all_steps)
               std::cout << "Local discrete maximum principle was satisfied at all time steps" << std::endl;
            else
               std::cout << "Local discrete maximum principle was NOT satisfied at all time steps" << std::endl;
         }
      }

      // evaluate errors for use in adaptive mesh refinement
      if (has_exact_solution) evaluate_error(cycle);
   }

   // output grid, solution, and viscosity and print convergence results
   output_results();
}

/** \brief 
 */
template<int dim>
void TransportProblem<dim>::transient_step(const Vector<double>       &old_solution_stage,
                                           const Vector<double>       &eval_solution_stage,
                                           const SparseMatrix<double> &mass_matrix,
                                           const SparseMatrix<double> &ss_matrix,
                                           const double               &dt,
                                           const double               &t_eval,
                                           const LinearSolver         &linear_solver)
{
   Vector<double> tmp_vector(n_dofs);

   // recompute the source if it is time-dependent
   if (source_time_dependent)
      assemble_ss_rhs(t_eval);

   // form transient rhs: system_rhs = M*u_old + dt*(ss_rhs - A*u_eval)
   system_rhs = 0;
   system_rhs.add(dt, ss_rhs);      //  now, system_rhs = dt*(ss_rhs)
   mass_matrix.vmult(tmp_vector, old_solution_stage);
   system_rhs.add(1.0, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs)
   ss_matrix.vmult(tmp_vector, eval_solution_stage);
   system_rhs.add(-dt, tmp_vector); //  now, system_rhs = M*u_old + dt*(ss_rhs - A*u_eval)
      
   // solve the linear system M*u_new = system_rhs
   system_matrix.copy_from(mass_matrix);
   apply_Dirichlet_BC(system_matrix, new_solution, system_rhs); // apply Dirichlet BC
   linear_solver.solve(system_matrix, system_rhs, new_solution);
     
   // distribute constraints
   constraints.distribute(new_solution);
}

/** \brief 
 */
template<int dim>
void TransportProblem<dim>::solve_steady_state()
{
   // entropy viscosity and FCT not yet implemented in steady-state
   Assert(parameters.scheme_option != 2, ExcNotImplemented());
   Assert(parameters.scheme_option != 3, ExcNotImplemented());

   // compute inviscid system matrix and steady-state right hand side (ss_rhs)
   assemble_inviscid_ss_matrix();
   assemble_ss_rhs(0.0);
   
   // add inviscid steady-state matrix to system matrix
   system_matrix.copy_from(inviscid_ss_matrix);

   if (parameters.scheme_option == 1)
   {
      // compute viscosity (both low-order and high-order)
      compute_low_order_viscosity();
      // compute viscous matrix
      compute_viscous_matrix(low_order_viscosity,low_order_viscous_matrix);
      // add low-order diffusion matrix to steady-state matrix
      system_matrix.add(1.0, low_order_viscous_matrix);
   }
   // enforce Dirichlet BC on total system matrix
   apply_Dirichlet_BC(system_matrix, new_solution, ss_rhs);
   // solve the linear system: ss_matrix*new_solution = ss_rhs
   //solve_linear_system(system_matrix, ss_rhs);
   LinearSolver linear_solver(parameters.linear_solver_option, constraints);
   linear_solver.solve(system_matrix, ss_rhs, new_solution);
   // distribute constraints
   constraints.distribute(new_solution);

   // compute bounds for maximum principle check
   compute_steady_state_max_principle_bounds();
   // check that solution satisfies maximum principle
   bool satisfied_max_principle = check_max_principle(0,0.0);
   // report if max principle was satisfied or not
   if (satisfied_max_principle)
      std::cout << "The local discrete maximum principle was satisfied." << std::endl;
   else
      std::cout << "The local discrete maximum principle was NOT satisfied." << std::endl;
}

/** \brief Output grid, solution, and viscosity to output file and print
 *         convergence table.
 */
template<int dim>
void TransportProblem<dim>::output_results()
{
   // output grid
   //------------
   if ((parameters.output_meshes) and (dim > 1))
      output_grid();

   // output solution
   output_solution(new_solution,
                   dof_handler,
                   "solution",
                   true);

   // write exact solution output file
   //---------------------------------
   if (parameters.output_exact_solution and has_exact_solution)
   {
      // create fine mesh on which to interpolate exact solution function
      Triangulation<dim> fine_triangulation;
      fine_triangulation.copy_triangulation(triangulation);
      unsigned int current_refinement_level = parameters.initial_refinement_level
         + parameters.n_refinement_cycles - 1;
      // define "fine" triangulation to be of refinement level 10, so refine
      // if refinement level is below this
      if (current_refinement_level < 10) fine_triangulation.refine_global(10-current_refinement_level);
      // create dof handler for fine mesh
      DoFHandler<dim> fine_dof_handler(fine_triangulation);
      fine_dof_handler.distribute_dofs(fe);
      // interpolate exact solution
      Vector<double> exact_solution_values(fine_dof_handler.n_dofs());
      exact_solution.set_time(parameters.end_time);
      VectorTools::interpolate(fine_dof_handler,
                               exact_solution,
                               exact_solution_values);

      // output exact solution to file
      output_solution(exact_solution_values,
                      fine_dof_handler,
                      "exact_solution",
                      false);
   }

   // write viscosity output file
   //----------------------------
   if (parameters.scheme_option != 0) {
      DataOut<dim> visc_out;
      visc_out.attach_dof_handler(dof_handler);
      // add viscosity data vector(s)
      switch (parameters.scheme_option)
      {
         case 1: {
            visc_out.add_data_vector(low_order_viscosity,"Low_Order_Viscosity",DataOut<dim>::type_cell_data);
            break;
         } case 2: {
            visc_out.add_data_vector(low_order_viscosity, "Low_Order_Viscosity", DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(entropy_viscosity,   "Entropy_Viscosity",   DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(high_order_viscosity,"High_Order_Viscosity",DataOut<dim>::type_cell_data);
            break;
         } case 3: {
            visc_out.add_data_vector(low_order_viscosity, "Low_Order_Viscosity", DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(entropy_viscosity,   "Entropy_Viscosity",   DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(high_order_viscosity,"High_Order_Viscosity",DataOut<dim>::type_cell_data);
            break;
         } default: {
            Assert(false, ExcNotImplemented());
            break;
         }
      }

      // determine output file extension
      std::string filename_extension;
      if (dim ==  1)
         filename_extension = ".gpl";
      else
         filename_extension = ".vtk";

      // create output filename
      std::stringstream viscosity_filename_ss;
      viscosity_filename_ss << "output/viscosity" << filename_extension;
      std::string viscosity_filename = viscosity_filename_ss.str();
      char *viscosity_filename_char = (char*)viscosity_filename.c_str();

      // create output filestream for exact solution
      std::ofstream viscosity_outstream(viscosity_filename_char);
      // build patches and write to file
      visc_out.build_patches(degree + 1);
      if (dim == 1)
         visc_out.write_gnuplot(viscosity_outstream);
      else
         visc_out.write_vtk(viscosity_outstream);
   }

   // output min and max bounds for DMP
   if (parameters.output_DMP_bounds and parameters.scheme_option == 3)
   {
      output_solution(min_values,dof_handler,"min_values",false);
      output_solution(max_values,dof_handler,"max_values",false);
   }

   // print convergence table
   //------------------------
   if (has_exact_solution) {
      convergence_table.set_precision("dx", 3);
      convergence_table.set_scientific("dx", true);
      if (not parameters.is_steady_state) {
         convergence_table.set_precision("dt", 3);
         convergence_table.set_scientific("dt", true);
         convergence_table.set_scientific("1/dt", true);
      }
      convergence_table.set_precision("L1 error", 3);
      convergence_table.set_scientific("L1 error", true);
      convergence_table.set_precision("L2 error", 3);
      convergence_table.set_scientific("L2 error", true);
      if (parameters.refinement_option == 2) {
         // evaluate temporal convergence rates
         convergence_table.evaluate_convergence_rates("L1 error", "1/dt", ConvergenceTable::reduction_rate_log2, 1);
         convergence_table.evaluate_convergence_rates("L2 error", "1/dt", ConvergenceTable::reduction_rate_log2, 1);
      } else {
         // evaluate spatial convergence rates
         convergence_table.evaluate_convergence_rates("L1 error", ConvergenceTable::reduction_rate_log2);
         convergence_table.evaluate_convergence_rates("L2 error", ConvergenceTable::reduction_rate_log2);
      }
      std::cout << std::endl;
      convergence_table.write_text(std::cout);
   }
}

/** \brief Outputs a solution to a file.
 *  \param [in] solution a vector of solution values.
 *  \param [in] dof_handler the dof handler associated with the solution values.
 *  \param [in] output_string string for the output filename.
 *  \param [in] append_viscosity the option to include a string for viscosity type in output filename
 */
template<int dim>
void TransportProblem<dim>::output_solution(const Vector<double>  &solution,
                                            const DoFHandler<dim> &dof_handler,
                                            const std::string     &output_string,
                                            const bool            &append_viscosity) const
{
   // create DataOut object for solution
   DataOut<dim> data_out;
   data_out.attach_dof_handler(dof_handler);
   data_out.add_data_vector(solution, "flux");
   data_out.build_patches(degree + 1);

   // create string for viscosity type if it is to be included in output filename
   std::string viscosity_string;
   if (append_viscosity)
   {
      switch (parameters.scheme_option) {
         case 0: {
            viscosity_string = "_galerkin";
            break;
         } case 1: {
            viscosity_string = "_low_order";
            break;
         } case 2: {
            viscosity_string = "_high_order";
            break;
         } case 3: {
            viscosity_string = "_FCT";
            break;
         } default: {
            Assert(false, ExcNotImplemented());
            break;
         }
      }
   } else {
      viscosity_string = "";
   }

   // create output filename for exact solution
   std::string filename_extension;
   if (dim == 1) filename_extension = ".gpl";
   else          filename_extension = ".vtk";

   std::stringstream filename_ss;
   filename_ss << "output/" << output_string << viscosity_string
      << "_" << parameters.problem_id << filename_extension;
   std::string filename = filename_ss.str();
   char *filename_char = (char*)filename.c_str();

   // create output filestream for exact solution
   std::ofstream output_filestream(filename_char);
   // write file
   if (dim == 1) data_out.write_gnuplot(output_filestream);
   else          data_out.write_vtk    (output_filestream);
}

/** \brief evaluate error between numerical and exact solution
 */
template<int dim>
void TransportProblem<dim>::evaluate_error(const unsigned int cycle)
{
   // assert that this function is only being called when an exact solution is available
   Assert(has_exact_solution,ExcInvalidState());

   // set time for exact solution function
   exact_solution.set_time(parameters.end_time);

   // error per cell
   Vector<double> difference_per_cell (n_cells);

   // compute L1 error
   VectorTools::integrate_difference(MappingQ<dim>(1),
                                     dof_handler,
                                     new_solution,
                                     exact_solution,
                                     difference_per_cell,
                                     cell_quadrature,
                                     VectorTools::L1_norm);
   const double L1_error = difference_per_cell.l1_norm();

   // compute L2 error
   VectorTools::integrate_difference(MappingQ<dim>(1),
                                     dof_handler,
                                     new_solution,
                                     exact_solution,
                                     difference_per_cell,
                                     cell_quadrature,
                                     VectorTools::L2_norm);
   const double L2_error = difference_per_cell.l2_norm();

   // compute average cell length
   const double avg_cell_length = std::pow((x_max-x_min),dim) / std::pow(n_cells,1.0/dim);

   // add error values to convergence table
   convergence_table.add_value("cycle", cycle);
   convergence_table.add_value("cells", n_cells);
   convergence_table.add_value("dofs", n_dofs);
   convergence_table.add_value("dx", avg_cell_length);
   if (not parameters.is_steady_state) {
      convergence_table.add_value("dt", dt_nominal);
      convergence_table.add_value("1/dt", 1.0/dt_nominal);
   }
   convergence_table.add_value("L1 error", L1_error);
   convergence_table.add_value("L2 error", L2_error);
}

/** \brief check that the local discrete max principle is satisfied.
 */
template<int dim>
bool TransportProblem<dim>::check_max_principle(const double &dt,
                                                const bool   &using_high_order)
{
   // machine precision for floating point comparisons
   const double machine_tolerance = 1.0e-12;

   // now set new precision
   std::cout.precision(15);

   // check that each dof value is bounded by its neighbors
   bool local_max_principle_satisfied = true;
   for (unsigned int i = 0; i < n_dofs; ++i) {
      // check if dof is a Dirichlet node or not - if it is, don't check max principle
      if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) == dirichlet_nodes.end())
      {
         double value_i = new_solution(i);
         // check lower bound
         if (value_i < min_values(i) - machine_tolerance) {
            local_max_principle_satisfied = false;
            // determine which condition was violated
            if (using_high_order)
            {
               std::cout << "      Max principle lower bound violated with dof "
                  << i << " of high-order solution: " << std::scientific
                  << value_i << " < " << min_values(i) << std::endl;
            } else {
               std::cout << "      Max principle lower bound violated with dof "
                  << i << " of low-order solution: " << std::scientific
                  << value_i << " < " << min_values(i) << std::endl;
               debug_max_principle_low_order (i,dt);
            }
         }
         // check upper bound
         if (value_i > max_values(i) + machine_tolerance) {
            local_max_principle_satisfied = false;
            // determine which condition was violated
            if (using_high_order)
            {
               std::cout << "      Max principle upper bound violated with dof "
                  << i << " of high-order solution: " << std::scientific
                  << value_i << " > " << max_values(i) << std::endl;
            } else {
               std::cout << "      Max principle upper bound violated with dof "
                  << i << " of low-order solution: " << std::scientific
                  << value_i << " > " << max_values(i) << std::endl;
               debug_max_principle_low_order (i,dt);
            }
         }
      }
   }

   // restore default precision and format
   std::cout.unsetf(std::ios_base::floatfield);

   // exit if max principle was violated
   if (not local_max_principle_satisfied)
   {
      std::cout << "Program terminated due to max-principle violation." << std::endl;
      std::exit(0);
   }

   return local_max_principle_satisfied;
}

/** \brief Debugging function used for determining why the maximum
 *         principle is failing for the low-order method;
 *         examines the conditions that are required to be met for
 *         the principle to be satisfied.
 */
template <int dim>
void TransportProblem<dim>::debug_max_principle_low_order(const unsigned int &i,
                                                          const double &dt)
{
   // flag to determine if no conditions were violated
   bool condition_violated = false;

   // get nonzero entries of row i of A
   std::vector<double>       row_values;
   std::vector<unsigned int> row_indices;
   unsigned int              n_col;
   get_matrix_row(low_order_ss_matrix,
                  i,
                  row_values,
                  row_indices,
                  n_col);

   // diagonal entry
   double Aii = low_order_ss_matrix(i,i);

   // check that system matrix diagonal elements are non-negative
   if (Aii < 0.0)
   {
      std::cout << "         DEBUG: diagonal element is negative: " << Aii << std::endl;
      condition_violated = true;
   }
   
   // loop over nonzero entries in row i to compute off-diagonal sum
   // for diagonal dominance check and also check that each off-diagonal
   // element is non-positive
   double off_diagonal_sum = 0.0;
   for (unsigned int k = 0; k < n_col; ++k)
   {
      unsigned int j = row_indices[k];
      double Aij = row_values[k];

      if (j != i)
      {
         // add to off-diagonal sum
         off_diagonal_sum += std::abs(Aij);
         // check that system matrix off-diagonal elements are non-positive
         if (Aij > 0.0)
         {
            std::cout << "         DEBUG: off-diagonal element (" << i << "," << j << ") is positive: " << Aij << std::endl;
            condition_violated = true;
         }
      }
   }

   // check that system matrix is diagonally dominant
   if (Aii < off_diagonal_sum)
   {
      std::cout << "         DEBUG: row is not diagonally dominant: Aii = "
         << Aii << ", off-diagonal sum = " << off_diagonal_sum << std::endl;
      condition_violated = true;
   }
     
   // check that CFL condition is satisfied if transient problem
   if (!parameters.is_steady_state)
   {
      double cfl = dt/lumped_mass_matrix(i,i)*low_order_ss_matrix(i,i);
      if (cfl > 1.0)
      {
         std::cout << "         DEBUG: row does not satisfy CFL condition: CFL = " << cfl << std::endl;
         condition_violated = true;
      }
   }

   // report if no conditions were violated
   if (not condition_violated)
      std::cout << "         DEBUG: No checks returned flags; deeper debugging is necessary." << std::endl;
}

/** \brief Computes min and max quantities for max principle
 */
template <int dim>
void TransportProblem<dim>::compute_max_principle_bounds(const double &dt)
{
   // initialize min and max values
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      max_values(i) = old_solution(i);
      min_values(i) = old_solution(i);
   }

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // find min and max values on cell
      double max_cell = old_solution(local_dof_indices[0]); // initialized to arbitrary value on cell
      double min_cell = old_solution(local_dof_indices[0]); // initialized to arbitrary value on cell
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
         double value_j = old_solution(local_dof_indices[j]);
         max_cell = std::max(max_cell, value_j);
         min_cell = std::min(min_cell, value_j);
      }

      // update the max and min values of neighborhood of each dof
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
         unsigned int i = local_dof_indices[j]; // global index
         max_values(i) = std::max(max_values(i), max_cell);
         min_values(i) = std::min(min_values(i), min_cell);
      }
   }

   // At this point, the min/max values of the old solution in the support
   // of test function i are stored in min_values(i) and max_values(i).
   // Now these values are multiplied by (1-dt/m(i))*sum_j(A(i,j)) and
   // added to dt/m(i)*b(i).

   // compute the upper and lower bounds for the maximum principle
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // compute sum of A_ij over row i
      // get nonzero entries of row i of A
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(low_order_ss_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);
      // add nonzero entries to get the row sum
      double row_sum = 0.0;
      for (unsigned int k = 0; k < n_col; ++k)
         row_sum += row_values[k];
      
      // compute the max and min values for the maximum principle
      max_values(i) = max_values(i)*(1.0 - dt/lumped_mass_matrix(i,i)*row_sum)
         + dt/lumped_mass_matrix(i,i)*ss_rhs(i);
      min_values(i) = min_values(i)*(1.0 - dt/lumped_mass_matrix(i,i)*row_sum)
         + dt/lumped_mass_matrix(i,i)*ss_rhs(i);
   }
}

/** \brief Computes min and max quantities for steady-state max principle
 */
template <int dim>
void TransportProblem<dim>::compute_steady_state_max_principle_bounds()
{
   // initialize min and max values
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      max_values(i) = new_solution(i);
      min_values(i) = new_solution(i);
   }

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // find min and max values on cell
      double max_cell = new_solution(local_dof_indices[0]); // initialized to arbitrary value on cell
      double min_cell = new_solution(local_dof_indices[0]); // initialized to arbitrary value on cell
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
         double value_j = new_solution(local_dof_indices[j]);
         max_cell = std::max(max_cell, value_j);
         min_cell = std::min(min_cell, value_j);
      }

      // update the max and min values of neighborhood of each dof
      for (unsigned int j = 0; j < dofs_per_cell; ++j) {
         unsigned int i = local_dof_indices[j]; // global index
         max_values(i) = std::max(max_values(i), max_cell);
         min_values(i) = std::min(min_values(i), min_cell);
      }
   }

   // compute the upper and lower bounds for the maximum principle
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // compute sum of A_ij over row i
      // get nonzero entries of row i of A
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(low_order_ss_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);
      // add nonzero entries to get the row sum
      double row_sum = 0.0;
      for (unsigned int k = 0; k < n_col; ++k)
         row_sum += row_values[k];
      
      // compute the max and min values for the maximum principle
      max_values(i) = max_values(i)*(1.0 - row_sum/low_order_ss_matrix(i,i))
         + ss_rhs(i)/low_order_ss_matrix(i,i);
      min_values(i) = min_values(i)*(1.0 - row_sum/low_order_ss_matrix(i,i))
         + ss_rhs(i)/low_order_ss_matrix(i,i);
   }
}

/** \brief Computes the high-order right hand side vector (G) and
 *         stores in system_rhs.
 */
template <int dim>
void TransportProblem<dim>::compute_high_order_rhs()
{
   // Form total system matrix, first adding the inviscid component
   high_order_ss_matrix.copy_from(inviscid_ss_matrix);
   // add viscous component to total system matrix (A)
   // Note that high_order_viscosity has already been computed in compute_low_order_viscosity(),
   // which was called for the assembly of the low-order system.
   compute_viscous_matrix(high_order_viscosity,high_order_viscous_matrix);
   // add viscous matrix to steady-state system matrix
   high_order_ss_matrix.add(1.0, high_order_viscous_matrix);

   // G is to be stored in system_rhs and is computed as G = ss_rhs - high_order_ss_matrix*old_solution
   system_rhs = ss_rhs;
   // now use ss_rhs as tmp vector to hold high_order_ss_matrix*old_solution
   high_order_ss_matrix.vmult(ss_rhs, old_solution);
   // complete computation of G
   system_rhs.add(-1.0, ss_rhs);
}

/** \brief Assembles the high-order coefficient matrix \f$\mathcal{A}\f$.
 *  \param [in] dt current time step size
 */
template <int dim>
void TransportProblem<dim>::assemble_flux_correction_matrix(const double &dt)
{
/*
   // reset A to zero
   flux_correction_matrix = 0;

   // add A1 matrix to A by looping over cells
   //--------------------------------------------------------
   std::vector<unsigned int> local_dof_indices(dofs_per_cell);
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (unsigned int i_cell = 0; cell != endc; ++cell, ++i_cell)
   {
      // get global dof indices for local dofs
      cell->get_dof_indices(local_dof_indices);

      // get low-order and high-order viscosities for cell
      double nu_L = low_order_viscosity (i_cell);
      double nu_H = high_order_viscosity(i_cell);
      // compute cell volume, needed for viscous bilinear form bK(phi_j,phi_i)
      double cell_volume = cell->measure();

      for (unsigned int ii = 0; ii < dofs_per_cell; ++ii)
      {
         // get global index for local dof i
         unsigned int i = local_dof_indices[ii];
         
         for (unsigned int jj = 0; jj < dofs_per_cell; ++jj)
         {
            // get global index for local dof j
            unsigned int j = local_dof_indices[jj];

            // compute viscous bilinear form bK(phi_j,phi_i)
            double bK;
            if (ii == jj)
               bK = cell_volume;
            else
               bK = cell_volume / (1.0 - dofs_per_cell);

            // compute cell contribution to A1(i,j)
            double A1_ij = dt*(nu_L - nu_H)*(old_solution(j) - old_solution(i))*bK;
            // add it to A(i,j)
            flux_correction_matrix.add(i,j,A1_ij);
         }
      }
   }
   
   // add A2 matrix to A by looping over all degrees of freedom
   //----------------------------------------------------------
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // Note that the high-order right hand side (G) is stored in system_rhs,
      // and auxiliary_mass_matrix is the "B" matrix

      // get nonzero entries in row i of B
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(auxiliary_mass_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);

      // loop over nonzero entries in row i of B
      for (unsigned int k = 0; k < n_col; ++k)
      {
         // get the column index of nonzero entry k
         unsigned int j = row_indices[k];

         // add A2(i,j) to A(i,j)
         flux_correction_matrix.add(i,j,dt*(auxiliary_mass_matrix(j,i)*system_rhs(i)
                                                  -auxiliary_mass_matrix(i,j)*system_rhs(j)));
      }
   }
*/
}

/** \brief Computes the limiting coefficient vectors \f$R^+\f$ and \f$R^-\f$,
 *         used for computing the high-order solution from the low-order
 *         solution.
 */
template <int dim>
void TransportProblem<dim>::compute_limiting_coefficients()
{
   // reset flux correction vector
   flux_correction_vector = 0;

   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // for Dirichlet nodes, set R_minus and R_plus to 1 because no limiting is needed
      if (std::find(dirichlet_nodes.begin(), dirichlet_nodes.end(), i) != dirichlet_nodes.end())
      {
         R_minus(i) = 1.0;
         R_plus(i) = 1.0;
      }
      else
      {
         // get nonzero entries in row i of high-order coefficient matrix A
         std::vector<double>       row_values;
         std::vector<unsigned int> row_indices;
         unsigned int              n_col;
         get_matrix_row(flux_correction_matrix,
                        i,
                        row_values,
                        row_indices,
                        n_col);
   
         double P_plus_i  = 0.0;
         double P_minus_i = 0.0;
         // compute P_plus, P_minus
         for (unsigned int k = 0; k < n_col; ++k)
         {
            // get value of nonzero entry k
            double Fij = row_values[k];
   
            P_plus_i  += std::max(0.0,Fij);
            P_minus_i += std::min(0.0,Fij);
         }
   
         // compute R_plus(i)
         if (P_plus_i != 0.0)
            R_plus(i) = std::min(1.0, Q_plus(i)/P_plus_i);
         else
            R_plus(i) = 1.0;
   
         // compute R_minus(i)
         if (P_minus_i != 0.0)
            R_minus(i) = std::min(1.0, Q_minus(i)/P_minus_i);
         else
            R_minus(i) = 1.0;
      }
   }

   // if user chose not to limit, set R+ and R- = 1 so that limiting
   // coefficients will equal 1 and thus no limiting will occur
   if (parameters.do_not_limit)
   {
      for (unsigned int i = 0; i < n_dofs; ++i)
      {
         R_minus(i) = 1.0;
         R_plus(i)  = 1.0;
      }
   }

   // compute limited flux correction sum
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // get values and indices of nonzero entries in row i of coefficient matrix A
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(flux_correction_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);

      // perform flux correction for dof i
      // Note that flux correction sum is sum_{j in I(S_i)} of Lij*Fij.
      for (unsigned int k = 0; k < n_col; ++k) {
         unsigned int j = row_indices[k];
         double Fij     = row_values[k];
         // compute limiting coefficient Lij
         double Lij;
         if (Fij >= 0.0)
            Lij = std::min(R_plus(i),R_minus(j));
         else
            Lij = std::min(R_minus(i),R_plus(j));

         // add Lij*Fij to flux correction sum
         flux_correction_vector(i) += Lij*Fij;
      }
   }
}

/** \brief Computes the high-order maximum-principle preserving solution
 *         using the low-order solution and the limiting coefficient
 *         matrix \f$\mathcal{L}\f$.
 */
/*
template <int dim>
void TransportProblem<dim>::compute_high_order_solution()
{
   
         // compute high-order solution for dof i
         new_solution(i) = new_solution(i) + flux_correction_sum / lumped_mass_matrix(i,i);
      }
      else // i corresponds to a Dirichlet node
      {
         // get Dirichlet value for node i
         std::map<unsigned int, double> boundary_values;
         VectorTools::interpolate_boundary_values(dof_handler,
                                            incoming_boundary,
                                            ConstantFunction<dim>(incoming_flux_value, 1),
                                            boundary_values);
         std::map<unsigned int, double>::iterator it = boundary_values.find(i);
         // ensure that Dirichlet node i was actually found in Dirichlet node list
         Assert(it != boundary_values.end(),ExcInvalidState());
         // set new solution to have Dirichlet value
         new_solution(i) = it->second;
      }
   }
}
*/

/** \brief Checks that the CFL condition is satisfied; If not, adjusts time step size.
    \param [in,out] dt time step size for current time step
 */
template <int dim>
double TransportProblem<dim>::enforce_CFL_condition(double &dt)
{
   // CFL is dt*speed/dx
   double max_speed_dx = 0.0; // max(speed/dx); max over all i of A(i,i)/mL(i,i)
   for (unsigned int i = 0; i < n_dofs; ++i)
      max_speed_dx = std::max(max_speed_dx, low_order_ss_matrix(i,i) / lumped_mass_matrix(i,i));

   // compute CFL number
   double proposed_CFL = dt * max_speed_dx;

   // if computed CFL number is greater than the set limit, then adjust dt
   double adjusted_CFL = proposed_CFL;
   if (proposed_CFL > parameters.CFL_limit)
   {
      std::cout << "CFL limit exceeded; dt has been adjusted" << std::endl;
      adjusted_CFL = parameters.CFL_limit;
      dt = adjusted_CFL / max_speed_dx;
   }

   // return adjusted CFL number
   return adjusted_CFL;
}

/** \brief Gets the values and indices of nonzero elements in a row of a sparse matrix.
 *  \param [in] matrix sparse matrix whose row will be retrieved
 *  \param [in] i index of row to be retrieved
 *  \param [out] row_values vector of values of nonzero entries of row i
 *  \param [out] row_indices vector of indices of nonzero entries of row i
 *  \param [out] n_col number of nonzero entries of row i
 */
template <int dim>
void TransportProblem<dim>::get_matrix_row(const SparseMatrix<double>      &matrix,
                                           const unsigned int              &i,
                                                 std::vector<double>       &row_values,
                                                 std::vector<unsigned int> &row_indices,
                                                 unsigned int              &n_col
                                          )
{
    // get first and one-past-last iterator for row
    SparseMatrix<double>::const_iterator matrix_iterator     = matrix.begin(i);
    SparseMatrix<double>::const_iterator matrix_iterator_end = matrix.end(i);

    // compute number of entries in row and then allocate memory
    n_col = matrix_iterator_end - matrix_iterator;
    row_values .resize(n_col);
    row_indices.resize(n_col);

    // loop over columns in row
    for(unsigned int k = 0; matrix_iterator != matrix_iterator_end; ++matrix_iterator, ++k)
    {
      row_values[k]  = matrix_iterator->value();
      row_indices[k] = matrix_iterator->column();
    }
}

/** \brief Gets a list of dofs subject to Dirichlet boundary conditions.
 *         This is necessary because max principle checks are not applied to these nodes.
 */
template <int dim>
void TransportProblem<dim>::get_dirichlet_nodes()
{
   // get map of Dirichlet dof indices to Dirichlet values
   std::map<unsigned int, double> boundary_values;
   VectorTools::interpolate_boundary_values(dof_handler,
                                            incoming_boundary,
                                            ConstantFunction<dim>(incoming_flux_value, 1),
                                            boundary_values);
   // extract dof indices from map
   dirichlet_nodes.clear();
   for (std::map<unsigned int, double>::iterator it = boundary_values.begin(); it != boundary_values.end(); ++it)
      dirichlet_nodes.push_back(it->first);
}
