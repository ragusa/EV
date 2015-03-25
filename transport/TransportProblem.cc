/** \brief constructor for TransportProblem class
 */
template<int dim>
TransportProblem<dim>::TransportProblem(const TransportParameters<dim> &parameters) :
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
      has_exact_solution(false),
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
   bool time_dependent = !parameters.is_steady_state;
   if (time_dependent)
      variables += ",t";

   // create function parser constants
   function_parser_constants["pi"]       = numbers::PI;
   function_parser_constants["x_min"]    = x_min;
   function_parser_constants["x_mid"]    = x_min + 0.5*(x_max-x_min);
   function_parser_constants["x_max"]    = x_max;

   // initialize exact solution function
   if (has_exact_solution)
      exact_solution_function.initialize(variables,
                                         exact_solution_string,
                                         function_parser_constants,
                                         time_dependent);

   // initialize initial conditions function
   if (!parameters.is_steady_state)
      initial_conditions.initialize(variables,
                                    initial_conditions_string,
                                    function_parser_constants,
                                    time_dependent);

   // initialize source function
   source_function.initialize(variables,
                              source_string,
                              function_parser_constants,
                              time_dependent);

   // initialize cross section function
   cross_section_function.initialize(variables,
                                     cross_section_string,
                                     function_parser_constants,
                                     time_dependent);

   // initialize Dirichlet boundary value function
   incoming_function.initialize(variables,
                                incoming_string,
                                function_parser_constants,
                                time_dependent);

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
      case 1: { // pure absorber
         Assert(dim < 3,ExcNotImplemented());

         x_min = 0.0;
         x_max = 1.0;

         incoming_string = "1";
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
            has_exact_solution = true;
            exact_solution_string = "if(x<=t,source/sigma + (incoming - source/sigma)*exp(-sigma*(x-x_min)),0)";
         }

         initial_conditions_string = "0";

         break;
      } case 2: { // void to absorber
         Assert(dim < 3,ExcNotImplemented());

         x_min = 0.0;
         x_max = 1.0;

         incoming_string = "1";
         function_parser_constants["incoming"]  = 1.0;

         cross_section_string = "if(x<x_mid, 0, 10)";
         function_parser_constants["sigma"]  = 10.0;

         source_time_dependent = false;
         source_string = "0";
         function_parser_constants["source"] = 0.0;

         has_exact_solution = true;
         if (parameters.is_steady_state) { // steady-state
            if (dim == 1)      // 1-D
               exact_solution_string = "if(x<x_mid, incoming, exp(-sigma*(x-x_mid)))";
            else if (dim == 2) // 2-D
               exact_solution_string = "if(y<x_mid, incoming, if(x<x_mid, incoming, exp(-sigma*(x-x_mid))))";
            else
               Assert(false,ExcNotImplemented());
         } else {                          // transient
            if (dim == 1)      // 1-D
               exact_solution_string = "if(x-t<x_min, if(x<x_mid, incoming, exp(-sigma*(x-x_mid))), 0)";
            else if (dim == 2) // 2-D
               exact_solution_string = "if(x-t<x_min, if(y<x_mid, incoming, if(x<x_mid, incoming, exp(-sigma*(x-x_mid)))), 0)";
            else
               Assert(false,ExcNotImplemented());
         }
         initial_conditions_string = "0";
         break;
      } case 3: { // void
         Assert(dim == 1,ExcNotImplemented());

         x_min = 0.0;
         x_max = 1.0;
         cross_section_string = "0";
         source_string = "0";
         source_time_dependent = false;

         incoming_string = "1";
         function_parser_constants["incoming"]  = 1.0;

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

         incoming_string = "0";

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

         incoming_string = "0";

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

         incoming_string = "0";

         cross_section_string = "1";
         has_exact_solution = true;
         exact_solution_string = "exp(-t)*sin(pi*x)"; // assume omega_x = 1 and c = 1
         source_string = "-exp(-t)*sin(pi*x) + pi*exp(-t)*cos(pi*x) + exp(-t)*sin(pi*x)";
         source_time_dependent = true;
         initial_conditions_string = exact_solution_string;
         break;
      } case 8: { // source in left half
         Assert(!parameters.is_steady_state,ExcNotImplemented());
         Assert(dim < 2,ExcNotImplemented());

         x_min = 0.0;
         x_max = 1.0;

         function_parser_constants["speed"] = 1.0;

         incoming_string = "0";
         function_parser_constants["incoming"] = 0.0;

         cross_section_string = "100";
         function_parser_constants["sigma"]  = 100.0;

         source_time_dependent = false;
         source_string = "if (x<x_mid,10,0)";
         function_parser_constants["source"] = 10.0;

         has_exact_solution = true;
         exact_solution_string = (std::string)"source/sigma*(1-exp(-sigma*" +
            "max(0,min(x,x_mid)-max(x-speed*t,0))))" +
            "*exp(-sigma*max(0,min(x,x_max)-max(x-speed*t,x_mid)))";

         initial_conditions_string = "0";
         break;
      } case 9: { // MMS-4
         Assert(dim == 1,ExcNotImplemented()); // assume 1-D
         Assert(not parameters.is_steady_state,ExcNotImplemented()); // assume not steady-state

         x_min = 0.0;
         x_max = 1.0;

         incoming_string = "0";

         cross_section_string = "1";
         has_exact_solution = true;
         exact_solution_string = "x*t"; // assume omega_x = 1 and c = 1
         source_string = "x + t + x*t";
         source_time_dependent = true;
         initial_conditions_string = exact_solution_string;
         break;
      } case 10: { // MMS-5
         Assert(dim == 1,ExcNotImplemented()); // assume 1-D

         x_min = 0.0;
         x_max = 1.0;

         if (parameters.is_steady_state) {
            incoming_string = "1";
            cross_section_string = "1";
            has_exact_solution = true;
            exact_solution_string = "1"; // assume omega_x = 1 and c = 1
            source_string = "1";
            source_time_dependent = false;
         } else {
            incoming_string = "t";
            cross_section_string = "1";
            has_exact_solution = true;
            exact_solution_string = "t"; // assume omega_x = 1 and c = 1
            source_string = "1 + t";
            source_time_dependent = true;
            initial_conditions_string = exact_solution_string;
         }

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
   low_order_diffusion_matrix .reinit(constrained_sparsity_pattern);
   consistent_mass_matrix   .reinit(constrained_sparsity_pattern);
   lumped_mass_matrix       .reinit(constrained_sparsity_pattern);
   if ((parameters.scheme_option == 2)||(parameters.scheme_option == 3)
      ||(parameters.scheme_option == 4))
      high_order_diffusion_matrix.reinit(constrained_sparsity_pattern);
      
   // reinitialize auxiliary matrices
   CompressedSparsityPattern compressed_unconstrained_sparsity_pattern(n_dofs);
   DoFTools::make_sparsity_pattern(dof_handler, compressed_unconstrained_sparsity_pattern);
   unconstrained_sparsity_pattern.copy_from(compressed_unconstrained_sparsity_pattern);
   viscous_bilinear_forms            .reinit(unconstrained_sparsity_pattern);
   // compute viscous bilinear forms
   compute_viscous_bilinear_forms();

   // reinitialize solution vector, system matrix, and rhs
   old_solution.reinit(n_dofs);
   older_solution.reinit(n_dofs);
   oldest_solution.reinit(n_dofs);
   old_stage_solution.reinit(n_dofs);
   new_solution.reinit(n_dofs);
   system_rhs  .reinit(n_dofs);
   ss_rhs      .reinit(n_dofs);

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
      
   // compute low-order viscosity
   compute_low_order_viscosity();
   // compute low-order diffusion matrix
   switch (parameters.low_order_diffusion_option) {
      case 1: {   // graph-Laplacian
         compute_graphLaplacian_diffusion_matrix(low_order_viscosity,low_order_diffusion_matrix);
         break;
      } case 2: { // standard Laplacian
         compute_standard_diffusion_matrix(low_order_viscosity,low_order_diffusion_matrix);
         break;
      } default: {
         Assert(false,ExcNotImplemented());
      }
   }
   // compute low-order steady-state matrix
   low_order_ss_matrix.copy_from(inviscid_ss_matrix);
   low_order_ss_matrix.add(1.0,low_order_diffusion_matrix);

   // enforce CFL condition on nominal dt size
   if (!parameters.is_steady_state)
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

/** \brief Computes a graph-Laplacian diffusion matrix.
 */
template <int dim>
void TransportProblem<dim>::compute_graphLaplacian_diffusion_matrix(
   const Vector<double> &viscosity,
   SparseMatrix<double> &diffusion_matrix)
{
   // reset viscous matrix
   diffusion_matrix = 0;

   // cell matrix
   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

   // FE values
   FEValues<dim> fe_values(fe, cell_quadrature, update_gradients | update_JxW_values);

   // loop over cells
   unsigned int i_cell = 0;
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell, ++i_cell) {
      // reset cell matrix to zero
      cell_matrix = 0;

      // reinitialize FE values
      fe_values.reinit(cell);

      // compute cell volume
      double cell_volume = cell->measure();

      // compute cell contribution to global system matrix
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
         for (unsigned int j = 0; j < dofs_per_cell; ++j) {
            double viscous_bilinear_form;
            if (j == i)
               viscous_bilinear_form = cell_volume;
            else
               viscous_bilinear_form = cell_volume/(1.0 - dofs_per_cell);
 
            cell_matrix(i,j) += viscosity(i_cell) * viscous_bilinear_form;
         }

      // aggregate local matrix and rhs to global matrix and rhs
      std::vector<unsigned int> local_dof_indices(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix,
                                             local_dof_indices,
                                             diffusion_matrix);

   }
}

/** \brief Computes a standard Laplacian diffusion matrix.
 */
template <int dim>
void TransportProblem<dim>::compute_standard_diffusion_matrix(
   const Vector<double> &viscosity,
   SparseMatrix<double> &diffusion_matrix)
{
   // reset viscous matrix
   diffusion_matrix = 0;

   // cell matrix
   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

   // FE values
   FEValues<dim> fe_values(fe, cell_quadrature, update_gradients | update_JxW_values);
   FEFaceValues<dim> fe_face_values(fe, face_quadrature, update_values |
      update_gradients | update_normal_vectors | update_JxW_values);

   // loop over cells
   unsigned int i_cell = 0;
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell, ++i_cell) {
      // reset cell matrix to zero
      cell_matrix = 0;

      // reinitialize FE values
      fe_values.reinit(cell);

      // compute cell contributions
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
               cell_matrix(i,j) += (
                  viscosity(i_cell) *
                  fe_values[flux].gradient(j, q) *
                  fe_values[flux].gradient(i, q)
                  ) * fe_values.JxW(q);

      // loop over faces of cell
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
         // if face is on outflow boundary
         if (cell->face(face)->at_boundary() && cell->face(face)->boundary_indicator() == 0) {
            // reinitialize FE face values
            fe_face_values.reinit(cell, face);

            for (unsigned int q = 0; q < n_q_points_face; ++q)
               for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                     cell_matrix(i,j) -= (
                        viscosity(i_cell) *
                        fe_face_values[flux].gradient(j, q) *
                        fe_face_values[flux].value(i, q) *
                        fe_face_values.normal_vector(q)
                     ) * fe_face_values.JxW(q);
         }

      // aggregate local matrix and rhs to global matrix and rhs
      std::vector<unsigned int> local_dof_indices(dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(cell_matrix,
                                             local_dof_indices,
                                             diffusion_matrix);

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
               cell->face(face)->set_boundary_indicator(1);
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

      switch (parameters.low_order_diffusion_option) {
         case 1: { // graph-Laplacian
            // initialize to zero for max() function
            double low_order_viscosity_cell = 0.0;
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
               for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  if (i != j)
                     low_order_viscosity_cell = std::max(low_order_viscosity_cell,
                        std::max(0.0,inviscid_ss_matrix(local_dof_indices[i],local_dof_indices[j]))/
                        (-viscous_bilinear_forms(local_dof_indices[i],local_dof_indices[j])));
      
            low_order_viscosity(i_cell) = low_order_viscosity_cell;
            break;
         } case 2: { // standard Laplacian
            low_order_viscosity(i_cell) = 0.5 * cell->measure();
            break;
         } default: {
            Assert(false,ExcNotImplemented());
         }
      }
   }
}

/** \brief run the problem
 */
template<int dim>
void TransportProblem<dim>::run()
{
   // initialize problem
   initialize_system();

   // create name of output subdirectory
   std::stringstream output_dir_ss;
   output_dir_ss << "output/problem_" << parameters.problem_id << "/";
   std::string output_dir = output_dir_ss.str();

   // determine time integrator string
   std::string timedisc_string;
   if (parameters.is_steady_state) timedisc_string = "ss";
   else
      switch (parameters.time_integrator_option) {
         case SSPRKTimeIntegrator<dim>::FE:
            {timedisc_string = "FE";      break;}
         case SSPRKTimeIntegrator<dim>::SSP2:
            {timedisc_string = "SSPRK22"; break;}
         case SSPRKTimeIntegrator<dim>::SSP3: 
            {timedisc_string = "SSPRK33"; break;}
         default: {ExcNotImplemented();}
      }

   // determine viscosity string
   std::string viscosity_string;
   switch (parameters.scheme_option) {
      case 0: {viscosity_string = "galerkin";   break;}
      case 1: {viscosity_string = "low_order";  break;}
      case 2: {viscosity_string = "high_order"; break;}
      case 3: {viscosity_string = "EVFCT";      break;}
      case 4: {viscosity_string = "GalFCT";     break;}
      default: {ExcNotImplemented();}
   }

   // create filename appendage
   std::stringstream appendage_ss;
   appendage_ss << "_" << parameters.problem_id << "_" << viscosity_string
      << "_" << timedisc_string;

   // create filename for exact solution
   std::stringstream filename_exact_ss;
   filename_exact_ss << "solution_" << parameters.problem_id << "_exact";

   // create post-processor object
   unsigned int final_refinement_level = parameters.initial_refinement_level
      + parameters.n_refinement_cycles - 1;
   PostProcessor<dim> postprocessor(parameters.output_meshes,
                                    parameters.output_exact_solution,
                                    parameters.save_convergence_results,
                                    has_exact_solution,
                                    exact_solution_function,
                                    parameters.end_time,
                                    dt_nominal,
                                    parameters.is_steady_state,
                                    parameters.refinement_mode,
                                    final_refinement_level,
                                    fe,
                                    output_dir,
                                    appendage_ss.str(),
                                    filename_exact_ss.str(),
                                    cell_quadrature);

   // create refinement handler object
   RefinementHandler<dim> refinement_handler(triangulation,
                                             n_cells,
                                             dt_nominal,
                                             parameters.refinement_mode,
                                             parameters.use_adaptive_refinement,
                                             parameters.time_refinement_factor);

   // loop over refinement cycles
   for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; ++cycle)
   {
      // refine
      refinement_handler.refine(cycle);

      // setup system - distribute finite elements, reintialize matrices and vectors
      setup_system();

      // update dt displayed in convergence table
      postprocessor.update_dt(dt_nominal);

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
                               domain_volume,
                               parameters.EV_time_discretization);

      // create linear solver object
      LinearSolver<dim> linear_solver(parameters.linear_solver_option,
                                      constraints,
                                      dof_handler,
                                      incoming_function);

      // create FCT object
      FCT<dim> fct(dof_handler,
                   lumped_mass_matrix,
                   consistent_mass_matrix,
                   linear_solver,
                   constrained_sparsity_pattern,
                   dirichlet_nodes,
                   n_dofs,
                   dofs_per_cell,
                   parameters.do_not_limit);

      // if problem is steady-state, then just do one solve; else loop over time
      if (parameters.is_steady_state) { // run steady-state problem
         solve_steady_state(linear_solver);
      } else {                          // run transient problem
         // interpolate initial conditions
         initial_conditions.set_time(0.0);
         VectorTools::interpolate(dof_handler,
                                  initial_conditions,
                                  new_solution);
         constraints.distribute(new_solution);

         // if last cycle, output initial conditions if user requested
         if ((cycle == parameters.n_refinement_cycles-1) and
            (parameters.output_initial_solution)) {
            std::stringstream IC_filename_ss;
            IC_filename_ss << "solution_" << parameters.problem_id << "_initial";
            postprocessor.output_solution(new_solution,
                                          dof_handler,
                                          IC_filename_ss.str());
         }

         // create SSP Runge-Kutta time integrator object
         SSPRKTimeIntegrator<dim> ssprk(parameters.time_integrator_option,
                                        n_dofs,
                                        linear_solver,
                                        constrained_sparsity_pattern);

         // time loop
         double t_new = 0.0;
         double t_old = 0.0;
         double old_dt   = dt_nominal;
         double older_dt = dt_nominal;
         old_solution    = new_solution;
         older_solution  = new_solution;
         oldest_solution = new_solution;
         const double t_end = parameters.end_time;
         bool in_transient = true;
         Vector<double> tmp_vector(n_dofs);
         unsigned int n = 1; // time step index
         while (in_transient)
         {
            // shorten dt if new time would overshoot end time, then update t_new
            double dt = dt_nominal;
            if (t_old + dt >= t_end) {
               dt = t_end - t_old;
               in_transient = false;
            }
            t_new = t_old + dt;
            std::cout << "   time step " << n << ": t = " << t_old << "->" << t_new << std::endl;

            // initialize SSPRK time step
            ssprk.initialize_time_step(old_solution,dt);

            switch (parameters.scheme_option)
            {
               case 0: { // unmodified Galerkin scheme

                  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
                  {
                     // get stage time
                     double t_stage = ssprk.get_stage_time();

                     // recompute steady-state rhs if it is time-dependent
                     if (source_time_dependent) {
                        assemble_ss_rhs(t_stage);
                     }

                     // advance by an SSPRK step
                     ssprk.step(consistent_mass_matrix,inviscid_ss_matrix,ss_rhs,true);
                  }
                  // retrieve the final solution
                  ssprk.get_new_solution(new_solution);

                  break;
               } case 1: { // solve low-order system

                  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
                  {
                     // get stage time
                     double t_stage = ssprk.get_stage_time();

                     // recompute steady-state rhs if it is time-dependent
                     if (source_time_dependent) {
                        assemble_ss_rhs(t_stage);
                     }

                     // advance by an SSPRK step
                     ssprk.step(lumped_mass_matrix,low_order_ss_matrix,ss_rhs,true);
                  }
                  // retrieve the final solution
                  ssprk.get_new_solution(new_solution);

                  break;
               } case 2: { // high-order system with entropy viscosity

                  // compute EV only at beginning of time step
                  if (parameters.EV_time_discretization != EntropyViscosity<dim>::FE) {
                     // recompute high-order steady-state matrix
                     entropy_viscosity = EV.compute_entropy_viscosity(
                        old_solution,
                        older_solution,
                        oldest_solution,
                        old_dt,
                        older_dt,
                        t_old);
                     for (unsigned int i_cell = 0; i_cell < n_cells; ++i_cell)
                        high_order_viscosity(i_cell) = std::min(
                           entropy_viscosity(i_cell),
                           low_order_viscosity(i_cell));
                     compute_graphLaplacian_diffusion_matrix(
                        high_order_viscosity,high_order_diffusion_matrix);
                     high_order_ss_matrix.copy_from(inviscid_ss_matrix);
                     high_order_ss_matrix.add(1.0,high_order_diffusion_matrix);
                  }

                  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
                  {
                     // get stage time
                     double t_stage = ssprk.get_stage_time();

                     // recompute steady-state rhs if it is time-dependent
                     if (source_time_dependent) {
                        assemble_ss_rhs(t_stage);
                     }

                     if (parameters.EV_time_discretization == EntropyViscosity<dim>::FE) {
                        // compute Galerkin solution
                        ssprk.step(consistent_mass_matrix,inviscid_ss_matrix,ss_rhs,false);
   
                        // get Galerkin solution
                        ssprk.get_intermediate_solution(new_solution);
   
                        // get old stage solution
                        ssprk.get_stage_solution(i,old_stage_solution);
   
                        // recompute high-order steady-state matrix
                        entropy_viscosity = EV.compute_entropy_viscosity(
                           new_solution,
                           old_stage_solution,
                           old_solution, // not used; supply arbitrary argument
                           dt,
                           old_dt,       // not used; supply arbitrary argument
                           t_stage);
                        for (unsigned int i_cell = 0; i_cell < n_cells; ++i_cell)
                           high_order_viscosity(i_cell) = std::min(
                              entropy_viscosity(i_cell),
                              low_order_viscosity(i_cell));
                        compute_graphLaplacian_diffusion_matrix(
                           high_order_viscosity,high_order_diffusion_matrix);
                        high_order_ss_matrix.copy_from(inviscid_ss_matrix);
                        high_order_ss_matrix.add(1.0,high_order_diffusion_matrix);
                     }

                     // advance by an SSPRK step
                     ssprk.step(consistent_mass_matrix,high_order_ss_matrix,ss_rhs,true);
                  }
                  // retrieve the final solution
                  ssprk.get_new_solution(new_solution);

                  break;
               } case 3: { // EV FCT
                  // assert that the graph-Laplacian low-order diffusion option was used
                  Assert(parameters.low_order_diffusion_option == 1, ExcNotImplemented());

                  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
                  {
                     // get stage time
                     double t_stage = ssprk.get_stage_time();

                     // recompute steady-state rhs if it is time-dependent
                     if (source_time_dependent) {
                        assemble_ss_rhs(t_stage);
                     }

                     // compute Galerkin solution
                     ssprk.step(consistent_mass_matrix,inviscid_ss_matrix,ss_rhs,false);

                     // get Galerkin solution
                     ssprk.get_intermediate_solution(new_solution);

                     // get old stage solution
                     ssprk.get_stage_solution(i,old_stage_solution);

                     // recompute high-order steady-state matrix
                     entropy_viscosity = EV.compute_entropy_viscosity(new_solution,old_stage_solution,old_solution,dt,old_dt,t_stage);
                     for (unsigned int i_cell = 0; i_cell < n_cells; ++i_cell)
                        high_order_viscosity(i_cell) = std::min(entropy_viscosity(i_cell),low_order_viscosity(i_cell));
                     compute_graphLaplacian_diffusion_matrix(high_order_viscosity,high_order_diffusion_matrix);
                     high_order_ss_matrix.copy_from(inviscid_ss_matrix);
                     high_order_ss_matrix.add(1.0,high_order_diffusion_matrix);

                     // advance by an SSPRK step
                     ssprk.step(consistent_mass_matrix,high_order_ss_matrix,ss_rhs,false);

                     // get old stage solution
                     ssprk.get_stage_solution(i,old_stage_solution);

                     // get intermediate solution
                     ssprk.get_intermediate_solution(new_solution);

                     // perform FCT
                     fct.solve_FCT_system(new_solution,
                                          old_stage_solution,
                                          low_order_ss_matrix,
                                          ss_rhs,
                                          dt,
                                          low_order_diffusion_matrix,
                                          high_order_diffusion_matrix);

                     // set stage solution to be FCT solution for this stage
                     ssprk.set_intermediate_solution(new_solution);

                     // finish computing stage solution
                     ssprk.complete_stage_solution();
                     
                  }
                  // retrieve the final solution
                  ssprk.get_new_solution(new_solution);

                  break;
               } case 4: { // Galerkin FCT
                  // assert that the graph-Laplacian low-order diffusion option was used
                  Assert(parameters.low_order_diffusion_option == 1, ExcNotImplemented());

                  for (unsigned int i = 0; i < ssprk.n_stages; ++i)
                  {
                     // get stage time
                     double t_stage = ssprk.get_stage_time();

                     // recompute steady-state rhs if it is time-dependent
                     if (source_time_dependent) {
                        assemble_ss_rhs(t_stage);
                     }

                     // compute Galerkin solution
                     ssprk.step(consistent_mass_matrix,inviscid_ss_matrix,ss_rhs,false);

                     // get Galerkin solution
                     ssprk.get_intermediate_solution(new_solution);

                     // get old stage solution
                     ssprk.get_stage_solution(i,old_stage_solution);

                     // perform FCT
                     fct.solve_FCT_system(new_solution,
                                          old_stage_solution,
                                          low_order_ss_matrix,
                                          ss_rhs,
                                          dt,
                                          low_order_diffusion_matrix,
                                          high_order_diffusion_matrix);

                     // set stage solution to be FCT solution for this stage
                     ssprk.set_intermediate_solution(new_solution);

                     // finish computing stage solution
                     ssprk.complete_stage_solution();
                     
                  }
                  // retrieve the final solution
                  ssprk.get_new_solution(new_solution);

                  break;
               } default: {
                  Assert(false, ExcNotImplemented());
                  break;
               }

            }

            // update old solution, time, and time step size
            oldest_solution = older_solution;
            older_solution  = old_solution;
            old_solution    = new_solution;
            older_dt = old_dt;
            old_dt   = dt;
            t_old = t_new;

            // increment time step index
            n++;
         }

         // report if DMP was satisfied at all time steps
         if ((parameters.scheme_option == 3)||(parameters.scheme_option == 4))
         {
            // report if local discrete maximum principle is satisfied at all time steps
            if (fct.check_DMP_satisfied())
               std::cout << "Local discrete maximum principle was satisfied at all time steps" << std::endl;
            else
               std::cout << "Local discrete maximum principle was NOT satisfied at all time steps" << std::endl;
         }
      }

      // evaluate errors for use in adaptive mesh refinement
      if (has_exact_solution) postprocessor.evaluate_error(new_solution,
                                                           dof_handler,
                                                           triangulation,
                                                           cycle);
   }

   // output grid and solution and print convergence results
   postprocessor.output_results(new_solution,dof_handler,triangulation);
   // output viscosities if they were used
   postprocessor.output_viscosity(low_order_viscosity,entropy_viscosity,high_order_viscosity,dof_handler);
}

/** \brief Solves the steady-state system.
 */
template<int dim>
void TransportProblem<dim>::solve_steady_state(LinearSolver<dim> &linear_solver)
{
   // entropy viscosity and FCT not yet implemented in steady-state
   Assert(parameters.scheme_option != 2, ExcNotImplemented());
   Assert(parameters.scheme_option != 3, ExcNotImplemented());
   Assert(parameters.scheme_option != 4, ExcNotImplemented());

   // compute inviscid system matrix and steady-state right hand side (ss_rhs)
   assemble_inviscid_ss_matrix();
   assemble_ss_rhs(0.0);
   
   // add inviscid steady-state matrix to system matrix
   system_matrix.copy_from(inviscid_ss_matrix);

   if (parameters.scheme_option == 1) // low-order scheme
   {
      // compute low-order viscosity
      compute_low_order_viscosity();
      // compute viscous matrix
      switch (parameters.low_order_diffusion_option) {
         case 1: {   // graph-Laplacian
            compute_graphLaplacian_diffusion_matrix(low_order_viscosity,low_order_diffusion_matrix);
            break;
         } case 2: { // standard Laplacian
            compute_standard_diffusion_matrix(low_order_viscosity,low_order_diffusion_matrix);
            break;
         } default: {
            Assert(false,ExcNotImplemented());
         }
      }
      // add low-order diffusion matrix to steady-state matrix
      system_matrix.add(1.0, low_order_diffusion_matrix);
   }
   // solve the linear system: ss_matrix*new_solution = ss_rhs
   linear_solver.solve(system_matrix, ss_rhs, new_solution);
   // distribute constraints
   constraints.distribute(new_solution);
}

/** \brief Checks that the CFL condition is satisfied; If not, adjusts time step size.
    \param [in,out] dt time step size for current time step
 */
template <int dim>
double TransportProblem<dim>::enforce_CFL_condition(double &dt)
{
   // CFL is dt*speed/dx
   double max_speed_dx = 0.0; // max(speed/dx); max over all i of A(i,i)/mL(i,i)
   for (unsigned int i = 0; i < n_dofs; ++i)
      max_speed_dx = std::max(max_speed_dx, std::abs(low_order_ss_matrix(i,i)) / lumped_mass_matrix(i,i));

   // compute CFL number
   double proposed_CFL = dt * max_speed_dx;

   // if computed CFL number is greater than the set limit, then adjust dt
   double adjusted_CFL = proposed_CFL;
   if (proposed_CFL > parameters.CFL_limit)
   {
      adjusted_CFL = parameters.CFL_limit;
      dt = adjusted_CFL / max_speed_dx;
   }

   // return adjusted CFL number
   return adjusted_CFL;
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
                                            1,
                                            ZeroFunction<dim>(),
                                            boundary_values);
   // extract dof indices from map
   dirichlet_nodes.clear();
   for (std::map<unsigned int, double>::iterator it = boundary_values.begin(); it != boundary_values.end(); ++it)
      dirichlet_nodes.push_back(it->first);
}
