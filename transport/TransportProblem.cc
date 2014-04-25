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
      cell_quadrature(degree+1),
      face_quadrature(degree+1),
      n_q_points_cell(cell_quadrature.size()),
      n_q_points_face(face_quadrature.size()),
      nonlinear_iteration(0),
      transport_direction(0.0),
      incoming_boundary(1),
      has_exact_solution(false),
      source_option(1),
      source_value(0.0),
      cross_section_option(1),
      cross_section_value(0.0),
      incoming_flux_value(0.0),
      is_linear(false),
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

   // determine if problem is linear or not; this will depend on which viscosity is used
   if (parameters.viscosity_option == 2) // entropy viscosity
      is_linear = false;
   else
      is_linear = true;

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
   std::map<std::string,double> function_parser_constants;
   function_parser_constants["sigma"] = cross_section_value;
   function_parser_constants["source"] = source_value;
   function_parser_constants["incoming"] = incoming_flux_value;
   function_parser_constants["pi"] = numbers::PI;

   // initialize exact solution function
   if (has_exact_solution)
   {
      exact_solution.initialize(variables,
                                exact_solution_string,
                                function_parser_constants,
                                true);
   }

   // initialize initial conditions function
   if (!parameters.is_steady_state)
   {
      initial_conditions.initialize(variables,
                                    initial_conditions_string,
                                    function_parser_constants,
                                    true);
   }
}

/** \brief process problem ID
 */
template<int dim>
void TransportProblem<dim>::process_problem_ID()
{
   switch (parameters.problem_id)
   {
      case 1: { // incoming flux, constant cross section, constant source
         cross_section_option = 1;
         cross_section_value = 0.3;
         source_option = 1;
         source_value = 1.2;
         incoming_flux_value = 3.0;
         has_exact_solution = true;
         Assert(dim < 3,ExcNotImplemented());
         exact_solution_string = "source/sigma + (incoming - source/sigma)*exp(-sigma*(x+1))";
         initial_conditions_string = exact_solution_string; // this creates pseudotransient
         break;
      } case 2: { // high cross section in half (1-D) or quarter (2-D) of domain, no source, incoming flux
         cross_section_option = 2;
         cross_section_value = 1.0e2;
         source_option = 1;
         source_value = 0.0e0;
         incoming_flux_value = 1.0e0;
         initial_conditions_string = "0";
         has_exact_solution = true;
         if (dim == 1)
            exact_solution_string = "if(x-t<-1, if(x<0, incoming, exp(-sigma*x)), 0)";
         else if (dim == 2)
            exact_solution_string = "if(x-t<-1, if(y<0, incoming, if(x<0.0, incoming, exp(-sigma*x))), 0)";
         else
            Assert(false,ExcNotImplemented());
         break;
      } case 3: {
         cross_section_option = 1;
         cross_section_value = 1.0e2;
         source_option = 2;
         source_value = 1.0e2;
         incoming_flux_value = 1.0e0;
         has_exact_solution = true;
         Assert(false,ExcNotImplemented());
         break;
      } case 4: {
         cross_section_option = 1;
         cross_section_value = 1.0e0;
         source_option = 3;
         source_value = 1.0e0;
         incoming_flux_value = 2.28735528717884;
         has_exact_solution = true;
         Assert(false,ExcNotImplemented());
         break;
      } case 5: { // linear advection problem with step IC
         cross_section_option = 1;
         cross_section_value = 0.0;
         source_option = 1;
         source_value = 0.0;
         incoming_flux_value = 1.0e0;
         has_exact_solution = true;
         if (parameters.is_steady_state)
            exact_solution_string = "1.0";
         else
            exact_solution_string = "if((x-t)<0,incoming,0)";
         initial_conditions_string = "if(x<0,incoming,0)";
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
   system_matrix         .reinit(constrained_sparsity_pattern);
   inviscid_system_matrix.reinit(constrained_sparsity_pattern);
   consistent_mass_matrix.reinit(constrained_sparsity_pattern);
   lumped_mass_matrix    .reinit(constrained_sparsity_pattern);
   if (parameters.viscosity_option == 4) // high-order max principle viscosity
   {
      auxiliary_mass_matrix        .reinit(constrained_sparsity_pattern);
      high_order_coefficient_matrix.reinit(constrained_sparsity_pattern);
   }
      

   // assemble mass matrices
   assemble_mass_matrices();

   // reinitialize auxiliary matrices
   CompressedSparsityPattern compressed_unconstrained_sparsity_pattern(n_dofs);
   DoFTools::make_sparsity_pattern(dof_handler, compressed_unconstrained_sparsity_pattern);
   unconstrained_sparsity_pattern.copy_from(compressed_unconstrained_sparsity_pattern);
   viscous_bilinear_forms            .reinit(unconstrained_sparsity_pattern);
   max_principle_viscosity_numerators.reinit(unconstrained_sparsity_pattern);
   // compute viscous bilinear forms
   compute_viscous_bilinear_forms();

   // reinitialize solution vector, system matrix, and rhs
   old_solution.reinit(n_dofs);
   new_solution.reinit(n_dofs);
   system_rhs  .reinit(n_dofs);
   ss_rhs      .reinit(n_dofs);
   R_plus      .reinit(n_dofs);
   R_minus     .reinit(n_dofs);

   // reinitialize max principle bounds
   min_values.reinit(n_dofs);
   max_values.reinit(n_dofs);

   // reinitialize viscosities
   entropy_viscosity        .reinit(triangulation.n_active_cells());
   low_order_viscosity      .reinit(triangulation.n_active_cells());
   high_order_viscosity     .reinit(triangulation.n_active_cells());
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

   // if using high-order max-principle viscosity, then assemble B-matrix
   //----------------------------
   if (parameters.viscosity_option == 4) // high-order max-principle viscosity
   {
      auxiliary_mass_matrix.copy_from(lumped_mass_matrix);
      auxiliary_mass_matrix.add(-1.0, consistent_mass_matrix);

      // at this point, B contains (M^L - M^C). Need to apply (M^L)^-1 to right.
      // This happens to be Bij = (M^L - M^C)ij / (M^L)jj, and this will
      // be achieved by iterating over the entries of B to divide by (M^L)jj
      SparseMatrix<double>::iterator matrix_iterator     = auxiliary_mass_matrix.begin();
      SparseMatrix<double>::iterator matrix_iterator_end = auxiliary_mass_matrix.end();
      for (; matrix_iterator != matrix_iterator_end; ++matrix_iterator)
      {
         unsigned int i = matrix_iterator->row();
         unsigned int j = matrix_iterator->column();
         double Bij = matrix_iterator->value() / lumped_mass_matrix(j,j);
         auxiliary_mass_matrix.set(i,j,Bij);
      }
   }
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

/** \brief Assemble the inviscid system matrix and steady-state right hand side.
 * 
 *         The "inviscid system matrix" will contain inviscid terms and Laplacian
 *         viscous terms if there are any (old, non-maximum-principle-preserving
 *         defininitions of viscosity).
 */
template<int dim>
void TransportProblem<dim>::assemble_system(const double &dt)
{
   inviscid_system_matrix = 0;
   ss_rhs = 0;
   max_principle_viscosity_numerators = 0;

   const TotalCrossSection<dim> total_cross_section(cross_section_option,cross_section_value);
   const TotalSource<dim>       total_source       (source_option,source_value);

   // FE values, for assembly terms
   FEValues<dim> fe_values(fe, cell_quadrature,
         update_values | update_gradients | update_quadrature_points
               | update_JxW_values);

   // FE face values, for boundary conditions
   FEFaceValues<dim> fe_face_values(fe, face_quadrature,
         update_values | update_quadrature_points | update_JxW_values
               | update_normal_vectors | update_gradients);

   FullMatrix<double> inviscid_cell_matrix(dofs_per_cell, dofs_per_cell);
   FullMatrix<double>          cell_matrix(dofs_per_cell, dofs_per_cell);
   Vector<double>              cell_rhs   (dofs_per_cell);

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // total cross section values at each quadrature point on cell
   std::vector<double> total_cross_section_values(n_q_points_cell);
   // total source values at each quadrature point on cell
   std::vector<double> total_source_values(n_q_points_cell);

   // get domain-averaged entropy if using entropy viscosity for this iteration
   max_entropy_deviation_domain = 0.0;
   bool need_to_compute_entropy_viscosity = (parameters.viscosity_option == 2) or
                                            (parameters.viscosity_option == 4);
   if (need_to_compute_entropy_viscosity)
      compute_entropy_domain_average();

   // cell iterator
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   // loop over cells
   unsigned int i_cell = 0;
   for (cell = dof_handler.begin_active(); cell != endc; ++cell, ++i_cell) {
      // initialize local matrix and rhs to zero
      inviscid_cell_matrix = 0;
      cell_matrix = 0;
      cell_rhs = 0;

      // reinitialize FE values
      fe_values.reinit(cell);

      // get total cross section for all quadrature points
      total_cross_section.value_list(fe_values.get_quadrature_points(),
                                     total_cross_section_values);
      // get total source for all quadrature points
      total_source.value_list(fe_values.get_quadrature_points(),
                              total_source_values);

      // compute viscosity
      double viscosity = compute_viscosity(cell, i_cell, fe_values, total_cross_section_values, total_source_values, dt);

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
               inviscid_cell_matrix(i,j) += (
                  // divergence term
                  fe_values[flux].value(i, q)
                     * transport_direction
                     * fe_values[flux].gradient(j, q) +
                  // total interaction term
                  fe_values[flux].value(i, q)
                     * total_cross_section_values[q]
                     * fe_values[flux].value(j, q)
                  ) * fe_values.JxW(q);

               // add to matrix
               cell_matrix(i, j) +=
                  // Laplacian viscosity term
                  viscosity
                     * fe_values[flux].gradient(i, q)
                     * fe_values[flux].gradient(j, q)
                     * fe_values.JxW(q);
            } // end j

            cell_rhs(i) +=
               // total source term
               fe_values[flux].value(i, q)
                  * total_source_values[q] * fe_values.JxW(q);
         } // end i
      } // end q

      // compute face contributions to global system;
      //  these arise due to integration by parts of Laplacian viscosity term
      // ------------------------------------------------------------------
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
      {
         if (cell->face(face)->at_boundary()) {
            fe_face_values.reinit(cell, face);

            for (unsigned int q = 0; q < n_q_points_face; ++q)
               for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                     cell_matrix(i, j) -= (viscosity
                           * fe_face_values.shape_value(i, q)
                           * fe_face_values.normal_vector(q)
                           * fe_face_values.shape_grad(j, q)
                           * fe_face_values.JxW(q));
         }
      } 

      // add inviscid terms to cell matrix
      cell_matrix.add(1.0, inviscid_cell_matrix);

      // aggregate local matrix and rhs to global matrix and rhs
      constraints.distribute_local_to_global(cell_matrix,
                                             cell_rhs,
                                             local_dof_indices,
                                             inviscid_system_matrix,
                                             ss_rhs);

      // aggregate local inviscid matrix into global matrix
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
         for (unsigned int j = 0; j < dofs_per_cell; ++j)
            max_principle_viscosity_numerators.add(local_dof_indices[i],
                                                   local_dof_indices[j],
                                                   inviscid_cell_matrix(i,j));
            
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

/** \brief computes the domain-averaged entropy and the max entropy
 *         deviation in the domain
 */
template <int dim>
void TransportProblem<dim>::compute_entropy_domain_average()
{
   FEValues<dim> fe_values(fe, cell_quadrature, update_values | update_JxW_values);

   // compute domain-averaged entropy
   //--------------------------------
   double domain_integral_entropy = 0.0;

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
      // reinitialize FE values
      fe_values.reinit(cell);
      // get solution values
      std::vector<double> new_solution_values(n_q_points_cell);
      fe_values[flux].get_function_values(new_solution,new_solution_values);
      // loop over quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         // compute entropy at quadrature point
         double entropy = 0.5 * new_solution_values[q] * new_solution_values[q];
         // add contribution of quadrature point to entropy integral
         domain_integral_entropy += entropy * fe_values.JxW(q);
      }
   }
   // domain-averaged entropy
   domain_averaged_entropy = domain_integral_entropy / domain_volume;

   // compute max deviation of entropy from domain-averaged entropy
   //--------------------------------------------------------------
   // loop over cells
   for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
      // reinitialize FE values
      fe_values.reinit(cell);
      // get old values and gradients
      std::vector<double> old_solution_values(n_q_points_cell);
      fe_values[flux].get_function_values(old_solution,old_solution_values);
      // loop over quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         // compute entropy at quadrature point
         double entropy = 0.5 * old_solution_values[q] * old_solution_values[q];
         // add contribution of quadrature point to entropy integral
         max_entropy_deviation_domain = std::max(max_entropy_deviation_domain,
               std::abs(entropy - domain_averaged_entropy));
      }
   }
}

/** \brief computes the viscosity in a cell
 */
template <int dim>
double TransportProblem<dim>::compute_viscosity(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                const unsigned int        &i_cell,
                                                const FEValues<dim>       &fe_values,
                                                const std::vector<double> &total_cross_section_values,
                                                const std::vector<double> &total_source_values,
                                                const double              &dt)
{
   double viscosity;
   // cell size
   double h = cell->diameter();
   // compute old definition of first-order viscosity
   double old_first_order_viscosity_cell = parameters.old_first_order_viscosity_coefficient * h;
   switch (parameters.viscosity_option) {
      case 0: { // no viscosity
         break;
      } case 1: { // old definition of first-order viscosity
         viscosity                   = old_first_order_viscosity_cell;
         low_order_viscosity(i_cell) = viscosity;
         break;
      } case 2: { // old definition of high-order viscosity
         low_order_viscosity(i_cell) = old_first_order_viscosity_cell;

         compute_entropy_viscosity(cell,
                                   i_cell,
                                   fe_values,
                                   total_cross_section_values,
                                   total_source_values,
                                   dt);

         // determine high-order viscosity: minimum of low-order viscosity and entropy viscosity
         viscosity = std::min(old_first_order_viscosity_cell, entropy_viscosity(i_cell));
         high_order_viscosity(i_cell) = viscosity;
         break;
      } case 3: { // maximum-principle preserving low-order viscosity
         // maximum-principle preserving viscosity does not add Laplacian term;
         // to avoid unnecessary branching in assembly loop, just add Laplacian
         // term with zero viscosity, so just keep viscosity with its initialization
         // of zero
         viscosity = 0.0;
         break;
      } case 4: { // maximum-principle preserving high-order viscosity
         compute_entropy_viscosity(cell,
                                   i_cell,
                                   fe_values,
                                   total_cross_section_values,
                                   total_source_values,
                                   dt);

         // maximum-principle preserving viscosity does not add Laplacian term;
         // to avoid unnecessary branching in assembly loop, just add Laplacian
         // term with zero viscosity, so just keep viscosity with its initialization
         // of zero
         viscosity = 0.0;
         break;
      } default: { // invalid viscosity option
         break;
      }
   }

   return viscosity;
}

/** \brief Computes the entropy viscosity in a cell.
 */
template <int dim>
void TransportProblem<dim>::compute_entropy_viscosity(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                      const unsigned int        &i_cell,
                                                      const FEValues<dim>       &fe_values,
                                                      const std::vector<double> &total_cross_section_values,
                                                      const std::vector<double> &total_source_values,
                                                      const double              &dt)
{
   // get previous time step (n) values and gradients
   std::vector<double>         new_values   (n_q_points_cell);
   std::vector<Tensor<1,dim> > new_gradients(n_q_points_cell);
   fe_values[flux].get_function_values   (new_solution, new_values);
   fe_values[flux].get_function_gradients(new_solution, new_gradients);
   // get previous previous time step (n-1) values
   std::vector<double>         old_values   (n_q_points_cell);
   fe_values[flux].get_function_values   (old_solution, old_values);

   // compute max entropy residual in cell
   //----------------------------------------------------------------------------
   // compute entropy values at each quadrature point on cell. The entropy definition s = 0.5*u^2 is used.
   std::vector<double> new_entropy_values(n_q_points_cell,0.0);
   std::vector<double> old_entropy_values(n_q_points_cell,0.0);
   for (unsigned int q = 0; q < n_q_points_cell; ++q) {
      new_entropy_values[q] = 0.5 * new_values[q] * new_values[q];
      old_entropy_values[q] = 0.5 * old_values[q] * old_values[q];
   }
   // compute entropy residual values at each quadrature point on cell
   std::vector<double> entropy_residual_values(n_q_points_cell,0.0);
   for (unsigned int q = 0; q < n_q_points_cell; ++q)
/*
      entropy_residual_values[q] = (new_entropy_values[q] - old_entropy_values[q])/dt
         + transport_direction * new_values[q] * new_gradients[q]
         + total_cross_section_values[q] * new_entropy_values[q];
*/
      entropy_residual_values[q] = (new_entropy_values[q] - old_entropy_values[q])/dt
         + transport_direction * new_values[q] * new_gradients[q]
         + total_cross_section_values[q] * new_values[q] * new_values[q]
         - total_source_values[q] * new_values[q];
   // determine maximum entropy residual in cell
   double max_entropy_residual = 0.0;
   for (unsigned int q = 0; q < n_q_points_cell; ++q) {
      max_entropy_residual  = std::max(max_entropy_residual,  std::abs(entropy_residual_values[q]));
   }

   // compute max jump in cell
   //----------------------------------------------------------------------------
   FEFaceValues<dim> fe_values_face         (fe, face_quadrature,
      update_values | update_gradients | update_JxW_values | update_normal_vectors);
   FEFaceValues<dim> fe_values_face_neighbor(fe, face_quadrature,
      update_values | update_gradients | update_JxW_values | update_normal_vectors);

   std::vector<double>         values_face            (n_q_points_face);
   std::vector<double>         values_face_neighbor   (n_q_points_face);
   std::vector<Tensor<1,dim> > gradients_face         (n_q_points_face);
   std::vector<Tensor<1,dim> > gradients_face_neighbor(n_q_points_face);
   std::vector<Point<dim> >    normal_vectors         (n_q_points_face);
   Vector<double> entropy_face(n_q_points_face);

   double max_jump_in_cell = 0.0;
   double max_jump_on_face = 0.0;
   for (unsigned int iface = 0; iface < faces_per_cell; ++iface)
   {
      typename DoFHandler<dim>::face_iterator face = cell->face(iface);
      if (face->at_boundary() == false)
      {
         Assert(cell->neighbor(iface).state() == IteratorState::valid, ExcInternalError());
         typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(iface);
         const unsigned int ineighbor = cell->neighbor_of_neighbor(iface);
         Assert(ineighbor < faces_per_cell, ExcInternalError());

         fe_values_face.reinit(         cell,    iface);
         fe_values_face_neighbor.reinit(neighbor,ineighbor);

         // get values on face
         fe_values_face.get_function_values(new_solution, values_face);

         // compute entropy at each quadrature point on face
         for (unsigned int q = 0; q < n_q_points_face; ++q)
            entropy_face[q] = 0.5 * values_face[q] * values_face[q];

         // get gradients on adjacent faces of current cell and neighboring cell
         fe_values_face.get_function_gradients(         new_solution, gradients_face);
         fe_values_face_neighbor.get_function_gradients(new_solution, gradients_face_neighbor);

         // get normal vectors
         normal_vectors = fe_values_face.get_normal_vectors();

         max_jump_on_face = 0.0;
         for (unsigned int q = 0; q < n_q_points_face; ++q)
         {
            double jump_dEdn = (gradients_face[q]*values_face[q]
               - gradients_face_neighbor[q]*values_face_neighbor[q]) * normal_vectors[q];
            double jump_on_face = std::abs(transport_direction * normal_vectors[q] * jump_dEdn);
            max_jump_on_face = std::max(max_jump_on_face, jump_on_face);
         }
      } // end if (at_boundary())
      max_jump_in_cell = std::max(max_jump_in_cell, max_jump_on_face);
   } // end face loop

   // compute entropy viscosity in cell
   //----------------------------------------------------------------------------
   entropy_viscosity(i_cell) = (parameters.entropy_viscosity_coefficient * max_entropy_residual
      + parameters.jump_coefficient * max_jump_in_cell) / max_entropy_deviation_domain;
}

/** \brief adds the viscous bilinear form for the maximum-principle preserving viscosity
 */
template <int dim>
void TransportProblem<dim>::add_viscous_matrix(const Vector<double> &viscosity)
{
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
                                             system_matrix);

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

/** \brief Computes the maximum-principle preserving low-order and high-order viscosities for each cell.
 */
template <int dim>
void TransportProblem<dim>::compute_max_principle_viscosity()
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
                  std::abs(max_principle_viscosity_numerators(local_dof_indices[i],local_dof_indices[j]))/
                  (-viscous_bilinear_forms(local_dof_indices[i],local_dof_indices[j])));
            }
         } 
      }
      low_order_viscosity(i_cell) = low_order_viscosity_cell;

      // compute the high-order viscosity
      high_order_viscosity(i_cell) = std::min(low_order_viscosity_cell, entropy_viscosity(i_cell));
   }
}

/** \brief solve the linear system
 */
template<int dim>
void TransportProblem<dim>::solve_linear_system(const SparseMatrix<double> &A,
                                                const Vector<double>       &b)
{
   switch (parameters.solver_option) {
      case 1: {
         SparseDirectUMFPACK A_direct;
         A_direct.initialize(A);
         A_direct.vmult(new_solution, b);
         break;
      }
      case 2: {
         SolverControl solver_control(1000, 1e-6);
         SolverBicgstab<> solver(solver_control);

         switch (parameters.preconditioner_option) {
            case 1: {
               solver.solve(A, new_solution, b,
                     PreconditionIdentity());
               break;
            }
            case 2: {
               PreconditionJacobi<> preconditioner;
               preconditioner.initialize(A, 1.0);
               solver.solve(A, new_solution, b,
                     preconditioner);
               break;
            }
            default: {
               Assert(false,ExcNotImplemented());
               break;
            }
         }
         break;
      }
      default: {
         Assert(false,ExcNotImplemented());
         break;
      }
   }

   constraints.distribute(new_solution);
}

/** \brief refine the grid
 */
template<int dim>
void TransportProblem<dim>::refine_grid() {
   if (parameters.use_adaptive_mesh_refinement) {
      Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

      KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim - 1>(3),
            typename FunctionMap<dim>::type(), new_solution,
            estimated_error_per_cell);

      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
            estimated_error_per_cell, 0.3, 0.03);

      triangulation.execute_coarsening_and_refinement();
   } else
      triangulation.refine_global(1);
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
      // generate mesh if in first cycle, else refine
      if (cycle == 0) {
         // domain is the dim-dimensional hypercube (-1,1)^dim
         GridGenerator::hyper_cube(triangulation, -1, 1);
         domain_volume = std::pow(2.0,dim);
         triangulation.refine_global(parameters.initial_refinement_level);
      } else {
         refine_grid();
      }

      // setup system - distribute finite elements, reintialize matrices and vectors
      setup_system();

      // print information
      std::cout << std::endl << "Cycle " << cycle << ':' << std::endl;
      std::cout << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl;
      std::cout << "   Number of degrees of freedom: " << n_dofs           << std::endl;

      // if problem is steady-state, then just do one solve; else loop over time
      if (parameters.is_steady_state) {
         // solve for steady-state solution
         solve_steady_state();
      } else {
         // have not yet added the capability for nonlinear transients
         Assert(is_linear,ExcNotImplemented());

         // interpolate initial conditions
         VectorTools::interpolate(dof_handler,
                                  initial_conditions,
                                  new_solution);
         // distribute hanging node and Dirichlet constraints to intial solution
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
         double t = 0.0;
         const double t_end = parameters.end_time;
         const double dt_const = parameters.time_step_size;
         const double machine_precision = 1.0e-15;
         bool in_transient = true;
         bool max_principle_satisfied = true;
         double old_dt = parameters.time_step_size; // time step size of previous time step, needed for entropy viscosity
         Vector<double> tmp_vector(n_dofs);
         unsigned int n = 1; // time step index
         while (in_transient)
         {
            // compute inviscid system matrix and steady-state right hand side (ss_rhs)
            // This inviscid system matrix will contain inviscid terms and Laplacian viscous terms if there are any.
            // Max-principle viscosity terms will be added in a separate step afterward.
            assemble_system(old_dt);
            // add inviscid component to total system matrix (A)
            system_matrix.copy_from(inviscid_system_matrix);
            // add viscous bilinear form for maximum-principle preserving viscosity
            bool using_max_principle_viscosity = (parameters.viscosity_option == 3) or (parameters.viscosity_option == 4);
            if (using_max_principle_viscosity) {
               // compute max-principle-preserving viscosity (both low-order and high-order)
               compute_max_principle_viscosity();
               // add viscous component to total system matrix (A)
               add_viscous_matrix(low_order_viscosity);
            }

            // update old solution to previous step's new solution
            // Note that this is done here because the old old_solution is needed in the
            // computation of the entropy viscosity time derivative term, which occurs in
            // assemble_system(), so the following update to old_solution needs to come
            // sometime after assemble_system() but before old_solution is used in the
            // assembly of the transient residual.
            old_solution = new_solution;

            // determine time step size
            double dt = parameters.time_step_size;
            if (t + dt_const > t_end) dt = t_end - t; // correct dt if new t would overshoot t_end 
            // enforce CFL condition and get CFL number
            double CFL;
            enforce_CFL_condition(dt,CFL);
            // increment time
            std::cout << "time step " << n << ": t = " << t << "->";
            t += dt;
            std::cout << t << ", CFL = " << CFL << std::endl;
            // determine if end of transient has been reached (within machine precision)
            in_transient = t < (t_end - machine_precision);

            // form transient rhs: system_rhs = M*u_old + dt*(ss_rhs - A*u_old)
            system_rhs = 0;
            system_rhs.add(dt, ss_rhs); //................. now, system_rhs = dt*(ss_rhs)
            lumped_mass_matrix.vmult(tmp_vector, old_solution);
            system_rhs.add(1.0, tmp_vector); //............ now, system_rhs = M*u_old + dt*(ss_rhs)
            system_matrix.vmult(tmp_vector, old_solution);
            system_rhs.add(-dt, tmp_vector); //............ now, system_rhs = M*u_old + dt*(ss_rhs - A*u_old), so it is complete
            // enforce the Dirichlet BC after solution has been formed
            apply_Dirichlet_BC(lumped_mass_matrix, new_solution, system_rhs);
            // solve the linear system M*u_new = system_rhs
            solve_linear_system(lumped_mass_matrix, system_rhs);
            // distribute constraints
            constraints.distribute(new_solution);

            // compute max principle min and max values
            compute_max_principle_bounds(dt);
            // check that local discrete maximum principle is satisfied at all time steps
            bool max_principle_satisfied_low_order = check_max_principle(dt,false);
            max_principle_satisfied = max_principle_satisfied and max_principle_satisfied_low_order;

            bool using_high_order = parameters.viscosity_option == 4;
            if (using_high_order) // high-order max-principle preserving viscosity
            {
               // check conservation of low-order solution
               double conservation_sum = 0.0;
               for (unsigned int i = 0; i < n_dofs; ++i)
                  conservation_sum += (new_solution(i) - old_solution(i))*lumped_mass_matrix(i,i);
               if (std::abs(conservation_sum) > 1.0e-15)
                  std::cout << "Low-order conservation sum != 0; sum = " << conservation_sum << std::endl;

               // compute high-order right hand side (G) and store in system_rhs
               compute_high_order_rhs();
               // compute high-order coefficient matrix (A)
               assemble_high_order_coefficient_matrix(dt);
               // compute the limiting coefficient matrix (L)
               compute_limiting_coefficients();
               // compute the high-order solution using the low-order solution and L
               compute_high_order_solution();
               // distribute constraints
               constraints.distribute(new_solution);

               // check that local discrete maximum principle is satisfied at all time steps
               bool max_principle_satisfied_high_order = check_max_principle(dt,true);
               max_principle_satisfied = max_principle_satisfied and max_principle_satisfied_high_order;
            }


            // update old time step size
            old_dt = dt;

            // increment time step index
            n++;
         }
         // report if local discrete maximum principle is satisfied at all time steps
         if (max_principle_satisfied)
            std::cout << "Local discrete maximum principle was satisfied at all time steps" << std::endl;
         else
            std::cout << "Local discrete maximum principle was NOT satisfied at all time steps" << std::endl;
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
void TransportProblem<dim>::solve_steady_state()
{
   // compute inviscid system matrix and steady-state right hand side (ss_rhs)
   // This inviscid system matrix will contain inviscid terms and Laplacian viscous terms if there are any.
   assemble_system(1.0e15);
   // add inviscid component to total system matrix (A)
   system_matrix.copy_from(inviscid_system_matrix);
   // add viscous bilinear form for maximum-principle preserving viscosity
   bool using_max_principle_viscosity = (parameters.viscosity_option == 3) or (parameters.viscosity_option == 4);
   if (using_max_principle_viscosity) {
      // compute max-principle-preserving viscosity (both low-order and high-order)
      compute_max_principle_viscosity();
      // add viscous component to total system matrix (A)
      add_viscous_matrix(low_order_viscosity);
   }
   // enforce Dirichlet BC on total system matrix
   apply_Dirichlet_BC(system_matrix, new_solution, ss_rhs);
   // solve the linear system: system_matrix*new_solution = ss_rhs
   solve_linear_system(system_matrix, ss_rhs);
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
   if (parameters.viscosity_option != 0) {
      DataOut<dim> visc_out;
      visc_out.attach_dof_handler(dof_handler);
      // add viscosity data vector(s)
      switch (parameters.viscosity_option)
      {
         case 1: {
            visc_out.add_data_vector(low_order_viscosity,"Old_Low_Order_Viscosity",DataOut<dim>::type_cell_data);
            break;
         } case 2: {
            visc_out.add_data_vector(low_order_viscosity, "Old_Low_Order_Viscosity",DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(entropy_viscosity,   "Entropy_Viscosity",        DataOut<dim>::type_cell_data);
            visc_out.add_data_vector(high_order_viscosity,"Old_High_Order_Viscosity",     DataOut<dim>::type_cell_data);
            break;
         } case 3: {
            visc_out.add_data_vector(low_order_viscosity,"Low_Order_Viscosity",DataOut<dim>::type_cell_data);
            break;
         } case 4: {
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
      viscosity_filename_ss << "output/viscosity_" << parameters.problem_id << filename_extension;
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

   // output min and max values for maximum principle check
   if ((parameters.viscosity_option == 3) or (parameters.viscosity_option == 4))
   {
      output_solution(min_values,dof_handler,"min_values",false);
      output_solution(max_values,dof_handler,"max_values",false);
   }

   // print convergence table
   //------------------------
   if (has_exact_solution) {
      convergence_table.set_precision("cell size", 3);
      convergence_table.set_scientific("cell size", true);
      convergence_table.set_precision("L2 error", 3);
      convergence_table.set_scientific("L2 error", true);
      convergence_table.evaluate_convergence_rates("L2 error", ConvergenceTable::reduction_rate_log2);
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
      switch (parameters.viscosity_option) {
         case 0: {
            viscosity_string = "_none";
            break;
         } case 1: {
            viscosity_string = "_old_first_order";
            break;
         } case 2: {
            viscosity_string = "_entropy";
            break;
         } case 3: {
            viscosity_string = "_low_order";
            break;
         } case 4: {
            viscosity_string = "_high_order";
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

   // error per cell
   Vector<double> difference_per_cell (triangulation.n_active_cells());

   // compute error with analytic solution
   VectorTools::integrate_difference(MappingQ<dim>(1),
                                     dof_handler,
                                     new_solution,
                                     exact_solution,
                                     difference_per_cell,
                                     QGauss<dim>(degree+1),
                                     VectorTools::L2_norm);

   // compute L2 error of vector of cell errors
   const double L2_error = difference_per_cell.l2_norm();

   // compute average cell length
   const unsigned int n_active_cells = triangulation.n_active_cells();
   const double avg_cell_length = std::pow(2.0,dim) / std::pow(n_active_cells,1.0/dim);

   // add error values to convergence table
   convergence_table.add_value("cycle", cycle);
   convergence_table.add_value("cells", n_active_cells);
   convergence_table.add_value("dofs", n_dofs);
   convergence_table.add_value("cell size", avg_cell_length);
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
               std::cout << "Max principle lower bound violated with dof "
                  << i << " of high-order solution: " << std::scientific
                  << value_i << " < " << min_values(i) << std::endl;
               debug_max_principle_high_order(i,dt);
            } else {
               std::cout << "Max principle lower bound violated with dof "
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
               std::cout << "Max principle upper bound violated with dof "
                  << i << " of high-order solution: " << std::scientific
                  << value_i << " > " << max_values(i) << std::endl;
               debug_max_principle_high_order(i,dt);
            } else {
               std::cout << "Max principle upper bound violated with dof "
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
   get_matrix_row(system_matrix,
                  i,
                  row_values,
                  row_indices,
                  n_col);

   // diagonal entry
   double Aii = system_matrix(i,i);

   // check that system matrix diagonal elements are non-negative
   if (Aii < 0.0)
   {
      std::cout << "   DEBUG: diagonal element is negative: " << Aii << std::endl;
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
            std::cout << "   DEBUG: off-diagonal element (" << i << "," << j << ") is positive: " << Aij << std::endl;
            condition_violated = true;
         }
      }
   }

   // check that system matrix is diagonally dominant
   if (Aii < off_diagonal_sum)
   {
      std::cout << "   DEBUG: row is not diagonally dominant: Aii = "
         << Aii << ", off-diagonal sum = " << off_diagonal_sum << std::endl;
      condition_violated = true;
   }
     
   // check that CFL condition is satisfied if transient problem
   if (!parameters.is_steady_state)
   {
      double cfl = dt/lumped_mass_matrix(i,i)*system_matrix(i,i);
      if (cfl > 1.0)
      {
         std::cout << "   DEBUG: row does not satisfy CFL condition: CFL = " << cfl << std::endl;
         condition_violated = true;
      }
   }

   // report if no conditions were violated
   if (not condition_violated)
      std::cout << "   DEBUG: No checks returned flags; deeper debugging is necessary." << std::endl;
}

/** \brief Debugging function used for determining why the maximum
 *         principle is failing for the low-order method;
 *         examines the conditions that are required to be met for
 *         the principle to be satisfied.
 */
template <int dim>
void TransportProblem<dim>::debug_max_principle_high_order(const unsigned int &i,
                                                           const double &dt)
{
   // flag to determine if no conditions were violated
   bool condition_violated = false;

   // check that the high order coefficient matrix A is skew symmetric
   // get nonzero entries of row i of A
   std::vector<double>       row_values;
   std::vector<unsigned int> row_indices;
   unsigned int              n_col;
   get_matrix_row(high_order_coefficient_matrix,
                  i,
                  row_values,
                  row_indices,
                  n_col);
   const double small_number = 1.0e-15;
   for (unsigned int k = 0; k < n_col; ++k)
   {
      unsigned int j = row_indices[k];
      if (j != i)
      {
         double Aij = row_values[k];
         double Aji = high_order_coefficient_matrix(j,i);
         // check that Aij = -Aji
         if (std::abs(Aij + Aji) > small_number)
         {
            std::cout << "   DEBUG: High-order coefficient matrix is not skew symmetric: A("
               << i << "," << j << ") = " << Aij << " != A(" << j << "," << i << ") = " << Aji
               << std::endl;
            condition_violated = true;
         }
      }
   }

   // check conservation of u^n -> uH^{n+1}
   // It is sufficient to check that sum_i G_i = 0 and sum_i B_ij = 0 (see Corollary 3.5)
   // first check sum_i G_i = 0
   double G_sum = 0.0;
   for (unsigned int i = 0; i < n_dofs; ++i)
      G_sum += system_rhs(i);
   if (std::abs(G_sum) > small_number)
   {
      std::cout << "   DEBUG: u^n -> uH^{n+1} is not conservative: sum_i G_i != 0; instead, = "
         << G_sum << std::endl;
      condition_violated = true;
   }

   // report if no conditions were violated
   if (not condition_violated)
      std::cout << "   DEBUG: No checks returned flags; deeper debugging is necessary." << std::endl;
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

   // compute the upper and lower bounds for the maximum principle
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // compute sum of A_ij over row i
      // get nonzero entries of row i of A
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(system_matrix,
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
      get_matrix_row(system_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);
      // add nonzero entries to get the row sum
      double row_sum = 0.0;
      for (unsigned int k = 0; k < n_col; ++k)
         row_sum += row_values[k];
      
      // compute the max and min values for the maximum principle
      max_values(i) = max_values(i)*(1.0 - row_sum/system_matrix(i,i))
         + ss_rhs(i)/system_matrix(i,i);
      min_values(i) = min_values(i)*(1.0 - row_sum/system_matrix(i,i))
         + ss_rhs(i)/system_matrix(i,i);
   }
}

/** \brief Computes the high-order right hand side vector (G) and
 *         stores in system_matrix.
 */
template <int dim>
void TransportProblem<dim>::compute_high_order_rhs()
{
   // Form total system matrix, first adding the inviscid component
   system_matrix.copy_from(inviscid_system_matrix);
   // add viscous component to total system matrix (A)
   // Note that high_order_viscosity has already been computed in compute_max_principle_viscosity(),
   // which was called for the assembly of the low-order system.
   add_viscous_matrix(high_order_viscosity);

   // G is to be stored in system_rhs and is computed as G = ss_rhs - system_matrix*old_solution
   system_rhs = ss_rhs;
   // now use ss_rhs as tmp vector to hold system_matrix*old_solution
   system_matrix.vmult(ss_rhs, old_solution);
   // complete computation of G
   system_rhs.add(-1.0, ss_rhs);
}

/** \brief Assembles the high-order coefficient matrix \f$\mathcal{A}\f$.
 *  \param [in] dt current time step size
 */
template <int dim>
void TransportProblem<dim>::assemble_high_order_coefficient_matrix(const double &dt)
{
   // reset A to zero
   high_order_coefficient_matrix = 0;

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
            high_order_coefficient_matrix.add(i,j,A1_ij);
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
         high_order_coefficient_matrix.add(i,j,dt*(auxiliary_mass_matrix(j,i)*system_rhs(i)
                                                  -auxiliary_mass_matrix(i,j)*system_rhs(j)));
      }
   }
}

/** \brief Computes the limiting coefficient vectors \f$R^+\f$ and \f$R^-\f$,
 *         used for computing the high-order solution from the low-order
 *         solution.
 */
template <int dim>
void TransportProblem<dim>::compute_limiting_coefficients()
{
   // small number around machine precision, used for floating point comparisons
   const double small_number = 1.0e-15;

   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // get nonzero entries in row i of high-order coefficient matrix A
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(high_order_coefficient_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);

      double P_plus  = 0.0;
      double P_minus = 0.0;
      // compute P_plus, P_minus
      for (unsigned int k = 0; k < n_col; ++k)
      {
         // get value of nonzero entry k
         double Aij = row_values[k];

         P_plus  += std::max(0.0,Aij);
         P_minus += std::min(0.0,Aij);
      }
      // compute Q_plus and Q_minus
      double Q_plus  = lumped_mass_matrix(i,i)*(max_values(i) - new_solution(i));
      double Q_minus = lumped_mass_matrix(i,i)*(min_values(i) - new_solution(i));

      // compute R_plus(i)
      if (std::abs(P_plus) > small_number)
         R_plus(i) = std::min(1.0, Q_plus/P_plus);
      else
         R_plus(i) = 1.0;

      // compute R_minus(i)
      if (std::abs(P_minus) > small_number)
         R_minus(i) = std::min(1.0, Q_minus/P_minus);
      else
         R_minus(i) = 1.0;
   }
}

/** \brief Computes the high-order maximum-principle preserving solution
 *         using the low-order solution and the limiting coefficient
 *         matrix \f$\mathcal{L}\f$.
 */
template <int dim>
void TransportProblem<dim>::compute_high_order_solution()
{
   for (unsigned int i = 0; i < n_dofs; ++i)
   {
      // get values and indices of nonzero entries in row i of coefficient matrix A
      std::vector<double>       row_values;
      std::vector<unsigned int> row_indices;
      unsigned int              n_col;
      get_matrix_row(high_order_coefficient_matrix,
                     i,
                     row_values,
                     row_indices,
                     n_col);

      // perform flux correction for dof i
      // Note that flux_correction_sum is sum_{j in I(S_i)} of Lij*Aij.
      double flux_correction_sum = 0.0;
      for (unsigned int k = 0; k < n_col; ++k) {
         unsigned int j = row_indices[k];
         double Aij     = row_values[k];
         if (Aij >= 0.0)
            flux_correction_sum += std::min(R_plus(i),R_minus(j))*Aij;
         else
            flux_correction_sum += std::min(R_minus(i),R_plus(j))*Aij;
      }

      // compute high-order solution for dof i
      new_solution(i) = new_solution(i) + flux_correction_sum / lumped_mass_matrix(i,i);
   }
}

/** \brief Checks that the CFL condition is satisfied; If not, adjusts time step size.
    \param [in,out] dt time step size for current time step
    \param [out] CFL CFL number
 */
template <int dim>
void TransportProblem<dim>::enforce_CFL_condition(double &dt, double &CFL) const
{
   // CFL is computed as max over all i of dt*A(i,i)/mL(i,i)
   double max_A_over_m = 0.0; // max over all i of A(i,i)/mL(i,i)
   for (unsigned int i = 0; i < n_dofs; ++i)
      max_A_over_m = std::max(max_A_over_m, system_matrix(i,i) / lumped_mass_matrix(i,i));

   // compute CFL number
   CFL = dt * max_A_over_m;

   // if computed CFL number is greater than the set limit, then adjust dt
   if (CFL > parameters.CFL_limit)
   {
      CFL = parameters.CFL_limit;
      dt = CFL / max_A_over_m;
   }
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
    row_values .reserve(n_col);
    row_indices.reserve(n_col);

    // loop over columns in row
    for(; matrix_iterator != matrix_iterator_end; ++matrix_iterator)
    {
      row_values .push_back(matrix_iterator->value());
      row_indices.push_back(matrix_iterator->column());
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
