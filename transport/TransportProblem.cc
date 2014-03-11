/** \brief constructor for TransportProblem class
 */
template<int dim>
TransportProblem<dim>::TransportProblem(const TransportParameters &parameters) :
      parameters(parameters),
      degree(parameters.degree),
      dof_handler(triangulation),
      fe(FE_Q<dim>(degree), 1),
      dofs_per_cell(fe.dofs_per_cell),
      cell_quadrature_formula(degree+1),
      face_quadrature_formula(degree+1),
      n_q_points_cell(cell_quadrature_formula.size()),
      n_q_points_face(face_quadrature_formula.size()),
      nonlinear_iteration(0),
      transport_direction(0.0)
{
}

/** \brief destructor for TransportProblem class
 */
template<int dim>
TransportProblem<dim>::~TransportProblem() {
   dof_handler.clear();
}

/** \brief set up the problem before assembly of the linear system
 */
template<int dim>
void TransportProblem<dim>::setup_system()
{
   dof_handler.distribute_dofs(fe);

   // reinitialize viscosity vectors
   max_viscosity          .reinit(triangulation.n_active_cells());
   entropy_viscosity      .reinit(triangulation.n_active_cells());
   max_principle_viscosity.reinit(triangulation.n_active_cells());

   // clear constraint matrix and make hanging node constraints for new mesh
   constraints.clear();
   DoFTools::make_hanging_node_constraints(dof_handler, constraints);
   constraints.close();

   // reinitialize sparsity pattern of system matrix
   CompressedSparsityPattern compressed_constrained_sparsity_pattern(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern(dof_handler,
                                               compressed_constrained_sparsity_pattern,
                                               constraints,
                                               false);
   constrained_sparsity_pattern.copy_from(compressed_constrained_sparsity_pattern);

   // allocate matrices to be used with maximum-principle preserving viscosity
   // reinitialize sparsity pattern of auxiliary matrices
   CompressedSparsityPattern compressed_unconstrained_sparsity_pattern(dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern(dof_handler, compressed_unconstrained_sparsity_pattern);
   unconstrained_sparsity_pattern.copy_from(compressed_unconstrained_sparsity_pattern);

   // reinitialize auxiliary matrices with sparsity pattern
   viscous_bilinear_forms            .reinit(unconstrained_sparsity_pattern);
   max_principle_viscosity_numerators.reinit(unconstrained_sparsity_pattern);

   // compute viscous bilinear forms
   compute_viscous_bilinear_forms();

   // reinitialize solution vector, system matrix, and rhs
   system_matrix   .reinit(constrained_sparsity_pattern);
   present_solution.reinit(dof_handler.n_dofs());
   system_rhs      .reinit(dof_handler.n_dofs());
}

/** \brief Computes viscous bilinear forms, to be used in the computation of
 *         maximum-principle preserving first order viscosity.
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

/** \brief assemble the system matrix and right hand side
 */
template<int dim>
void TransportProblem<dim>::assemble_system()
{
   max_principle_viscosity_numerators = 0;

   const TotalCrossSection<dim> total_cross_section(parameters);
   const TotalSource<dim> total_source(parameters);

   // FE values, for assembly terms
   FEValues<dim> fe_values(fe, cell_quadrature_formula,
         update_values | update_gradients | update_quadrature_points
               | update_JxW_values);

   // FE face values, for boundary conditions
   FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
         update_values | update_quadrature_points | update_JxW_values
               | update_normal_vectors | update_gradients);

   FullMatrix<double> inviscid_cell_matrix(dofs_per_cell, dofs_per_cell);
   FullMatrix<double>          cell_matrix(dofs_per_cell, dofs_per_cell);
   Vector<double> cell_rhs(dofs_per_cell);

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // total cross section values for an energy group and direction at each
   //  quadrature point on cell
   std::vector<double> total_cross_section_values(n_q_points_cell);
   // total source values for an energy group and direction at each
   //  quadrature point on cell
   std::vector<double> total_source_values(n_q_points_cell);

   // create extractor
   FEValuesExtractors::Scalar flux(0);

   // cell iterator
   typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active(), endc = dof_handler.end();

   // compute domain volume for denominator of domain-averaged entropy
   double domain_volume = 0.0;
   for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
      // reinitialize FE values
      fe_values.reinit(cell);
      // loop over quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         domain_volume += fe_values.JxW(q);
   }

   // get domain-averaged entropy if using entropy viscosity for this iteration
   // ---------------------------------------------------------------------
   double domain_averaged_entropy;
   double max_entropy_deviation_domain = 0.0;
   if ((parameters.viscosity_type == 2)&&(nonlinear_iteration != 0)) {
      // compute domain-averaged entropy
      double domain_integral_entropy = 0.0;
      // loop over cells
      for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
         // reinitialize FE values
         fe_values.reinit(cell);
         // loop over quadrature points
         for (unsigned int q = 0; q < n_q_points_cell; ++q) {
            // get angular flux at quadrature point
            double angular_flux = 0.0;
            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
               angular_flux += fe_values[flux].value(j, q);
            }
            // compute entropy at quadrature point
            double entropy = 1.0/2.0 * angular_flux * angular_flux;
            // add contribution of quadrature point to entropy integral
            domain_integral_entropy += entropy * fe_values.JxW(q);
         }
      }
      // domain-averaged entropy
      domain_averaged_entropy = domain_integral_entropy / domain_volume;

      // find max deviation of entropy from domain entropy average
      // loop over cells
      for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
         // reinitialize FE values
         fe_values.reinit(cell);
         // get old values and gradients
         std::vector<double> old_values(n_q_points_cell);
         fe_values[flux].get_function_values(old_solution,old_values);
         // loop over quadrature points
         for (unsigned int q = 0; q < n_q_points_cell; ++q) {
            // get angular flux at quadrature point
            double angular_flux = 0.0;
            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
               angular_flux += old_values[q];
            }
            // compute entropy at quadrature point
            double entropy = 1.0/2.0 * angular_flux * angular_flux;
            // add contribution of quadrature point to entropy integral
            max_entropy_deviation_domain = std::max(max_entropy_deviation_domain,
                  std::abs(entropy-domain_averaged_entropy));
         }
      }
   }

   // loop over cells
   unsigned int i_cell = 0;
   for (cell = dof_handler.begin_active(); cell != endc; ++cell, ++i_cell) {
      // initialize local matrix and rhs to zero
      inviscid_cell_matrix = 0;
      cell_matrix = 0;
      cell_rhs = 0;

      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // reinitialize FE values
      fe_values.reinit(cell);

      // get total cross section for all quadrature points
      total_cross_section.value_list(fe_values.get_quadrature_points(),
                                     total_cross_section_values);
      // get total source for all quadrature points
      total_source.value_list(fe_values.get_quadrature_points(),
                              total_source_values);

      // compute viscosity
      // ------------------------------------------------------------------
      double viscosity;
      double h = cell->diameter();
      double max_viscosity_cell = parameters.max_viscosity_coefficient * h;
      max_viscosity(i_cell) = max_viscosity_cell;
      switch (parameters.viscosity_type) {
         case 0: {
            break;
         } case 1: {
            viscosity = max_viscosity_cell;
            break;
         } case 2: {
            if (nonlinear_iteration != 0) {
               // get old values and gradients
               std::vector<double> old_values(n_q_points_cell);
               std::vector<Tensor<1,dim> > old_gradients(n_q_points_cell);
               fe_values[flux].get_function_values(old_solution,old_values);
               fe_values[flux].get_function_gradients(old_solution,old_gradients);

               // compute entropy values at each quadrature point on cell
               std::vector<double> entropy_values(n_q_points_cell,0.0);
               for (unsigned int q = 0; q < n_q_points_cell; ++q)
                  entropy_values[q] = 1.0/2.0 * old_values[q] * old_values[q];
               // compute entropy residual values at each quadrature point on cell
               std::vector<double> entropy_residual_values(n_q_points_cell,0.0);
               for (unsigned int q = 0; q < n_q_points_cell; ++q)
                  entropy_residual_values[q] = std::abs(transport_direction *
                        old_values[q] * old_gradients[q]
                                                        + total_cross_section_values[q] * entropy_values[q]);
               // compute entropy deviation values at each quadrature point on cell
               std::vector<double> entropy_deviation_values(n_q_points_cell,0.0);
               for (unsigned int q = 0; q < n_q_points_cell; ++q)
                  entropy_deviation_values[q] = std::abs(entropy_values[q] - domain_averaged_entropy);
               // determine cell maximum entropy residual and entropy deviation
               double max_entropy_residual = 0.0;
               for (unsigned int q = 0; q < n_q_points_cell; ++q) {
                  max_entropy_residual = std::max(max_entropy_residual,entropy_residual_values[q]);
               }

               // compute entropy viscosity
               double entropy_viscosity_cell = parameters.entropy_viscosity_coefficient *
                     h * h * max_entropy_residual / max_entropy_deviation_domain;
               entropy_viscosity(i_cell) = entropy_viscosity_cell;

               // determine viscosity: minimum of first-order viscosity and entropy viscosity
               viscosity = std::min(max_viscosity_cell, entropy_viscosity_cell);
            } else
               viscosity = max_viscosity_cell;

            break;
         } case 3: {
            // maximum-principle preserving viscosity does not add Laplacian term;
            // to avoid unnecessary branching in assembly loop, just add Laplacian
            // term with zero viscosity, so just keep viscosity with its initialization
            // of zero
            viscosity = 0.0;
            break;
         } default: {
            Assert(false,ExcNotImplemented());
            break;
         }
      }
      // compute cell contributions to global system
      // ------------------------------------------------------------------
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
                  // viscosity term
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

      // add face terms; these arise due to integration by parts of viscosity term
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
                                             system_matrix,
                                             system_rhs);

      // aggregate local inviscid matrix global matrix
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
         for (unsigned int j = 0; j < dofs_per_cell; ++j)
            max_principle_viscosity_numerators.add(local_dof_indices[i],
                                                   local_dof_indices[j],
                                                   inviscid_cell_matrix(i,j));
            
   } // end cell

   // add viscous bilinear form for maximum-principle preserving viscosity
   // ---------------------------------------------------------------------------
   if (parameters.viscosity_type == 3) {
      // compute maximum-principle preserving viscosity
      compute_max_principle_viscosity();

      unsigned int i_cell = 0;
      for (cell = dof_handler.begin_active(); cell != endc; ++cell, ++i_cell) {
         // reset cell matrix to zero
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
                  viscous_bilinear_form = -1.0/(dofs_per_cell - 1.0)*cell_volume;
 
               cell_matrix(i,j) += max_principle_viscosity(i_cell) * viscous_bilinear_form;
            }
         }

         // aggregate local matrix and rhs to global matrix and rhs
         cell->get_dof_indices(local_dof_indices);
         constraints.distribute_local_to_global(cell_matrix,
                                                local_dof_indices,
                                                system_matrix);

      }
   }

   // apply boundary conditions
   // ---------------------------------------------------------------------------
   // reset boundary indicators to zero
   for (cell = dof_handler.begin_active(); cell != endc; ++cell)
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
            ++face)
         if (cell->face(face)->at_boundary())
            cell->face(face)->set_boundary_indicator(0);

   // loop over cells
   for (cell = dof_handler.begin_active(); cell != endc; ++cell) {
      // loop over faces of cell
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
            ++face) {
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
   // apply Dirichlet boundary condition
   std::map<unsigned int, double> boundary_values;
   VectorTools::interpolate_boundary_values(dof_handler,
                                            1,
                                            ConstantFunction<dim>(parameters.incoming_flux, 1),
                                            boundary_values);
   MatrixTools::apply_boundary_values(boundary_values,
                                      system_matrix,
                                      present_solution,
                                      system_rhs);

} // end assembly

/** \brief Computes the maximum-principle preserving first order viscosity for each cell.
 */
template <int dim>
void TransportProblem<dim>::compute_max_principle_viscosity()
{
   // local dof indices
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   // loop over cells to compute first order viscosity at each quadrature point
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   unsigned int i_cell = 0;
   for (; cell != endc; ++cell, ++i_cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         max_principle_viscosity(i_cell) = 0.0;
         for (unsigned int i = 0; i < dofs_per_cell; ++i) {
            for (unsigned int j = 0; j < dofs_per_cell; ++j) {
               if (i != j) {
                  max_principle_viscosity(i_cell) = std::max(max_principle_viscosity(i_cell),
                     std::abs(max_principle_viscosity_numerators(local_dof_indices[i],local_dof_indices[j]))/
                     (-viscous_bilinear_forms(local_dof_indices[i],local_dof_indices[j])));
               }
            }
         }
      }
   }
}

/** \brief solve the linear system
 */
template<int dim>
void TransportProblem<dim>::solve_linear_system() {
   switch (parameters.solver_option) {
      case 1: {
         SparseDirectUMFPACK A_direct;
         A_direct.initialize(system_matrix);
         A_direct.vmult(present_solution, system_rhs);
         break;
      }
      case 2: {
         SolverControl solver_control(1000, 1e-6);
         SolverBicgstab<> solver(solver_control);

         switch (parameters.preconditioner_option) {
            case 1: {
               solver.solve(system_matrix, present_solution, system_rhs,
                     PreconditionIdentity());
               break;
            }
            case 2: {
               PreconditionJacobi<> preconditioner;
               preconditioner.initialize(system_matrix, 1.0);
               solver.solve(system_matrix, present_solution, system_rhs,
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

   constraints.distribute(present_solution);
}

/** \brief refine the grid
 */
template<int dim>
void TransportProblem<dim>::refine_grid() {
   if (parameters.use_adaptive_mesh_refinement) {
      Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

      KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim - 1>(3),
            typename FunctionMap<dim>::type(), present_solution,
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
   // decide that transport direction is the unit x vector
   transport_direction[0] = 1.0;

   // loop over refinement cycles
   for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; ++cycle)
   {
      // generate mesh if in first cycle, else refine
      if (cycle == 0) {
         GridGenerator::hyper_cube(triangulation, -1, 1);
         triangulation.refine_global(parameters.initial_refinement_level);
      } else {
         refine_grid();
      }

      // setup system - distribute finite elements, reintialize matrices and vectors
      setup_system();

      // print information
      std::cout << std::endl << "Cycle " << cycle << ':' << std::endl;
      std::cout << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl;
      std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()           << std::endl;

      // if problem is steady-state, then just do one solve; else loop over time
      if (parameters.is_steady_state) {
         // solve for steady-state solution
         solve_step();
      } else {
         // time loop
         double t = 0.0;
         const double t_end = parameters.end_time;
         const double dt_const = parameters.time_step_size;
         const double machine_precision = 1.0e-15;
         bool in_transient = true;
         while (in_transient)
         {
            // determine time step size
            double dt = parameters.time_step_size;
            if (t + dt_const > t_end) dt = t_end - t; // correct dt if new t would overshoot t_end 
            // increment time
            t += dt;
            // determine if end of transient has been reached (within machine precision)
            in_transient = t < t_end - machine_precision;

            // solve for current time solution
            solve_step();
         }
      }

      // evaluate errors for use in adaptive mesh refinement
      if (parameters.exact_solution_id != 0) evaluate_error(cycle);
   }

   // output grid, solution, and viscosity and print convergence results
   output_results();
}

/** \brief 
 */
template<int dim>
void TransportProblem<dim>::solve_step()
{
   // solve system. if using entropy viscosity, begin nonlinear iteration, else
   // just solve linear system
   nonlinear_iteration = 0;
   if (parameters.viscosity_type == 2) {
      bool converged = false; // converged nonlinear iteration
      for (unsigned int iter = 0; iter < parameters.max_nonlinear_iterations; ++iter) {
         std::cout << "   Nonlinear iteration " << iter;
         if (iter == 0) std::cout << std::endl;
         assemble_system();
         solve_linear_system();
         // if not the first iteration, evaluate the convergence criteria
         if (nonlinear_iteration != 0) {
            // evaluate the difference between the current and previous solution iterate
            double old_norm = old_solution.l2_norm();
            old_solution -= present_solution;
            double difference_norm = old_solution.l2_norm();
            double relative_difference = difference_norm / old_norm;
            std::cout << ": Error: " << relative_difference << std::endl;
            if (relative_difference < parameters.relative_difference_tolerance) {
               converged = true;
               break;
            }
         }
         // update the old solution and iteration number
         old_solution = present_solution;
         nonlinear_iteration++;
      }
      // report if the solution did not converge
      if (!converged) {
         std::cout << "The solution did not converge in " << parameters.max_nonlinear_iterations << " iterations";
         std::cout << std::endl;
      }
   } else {
      // system is linear and requires just one solve
      assemble_system();
      solve_linear_system();
   }
   // check that solution is non-negative
   check_solution_nonnegative();
}

/** \brief Output grid, solution, and viscosity to output file and print
 *         convergence table.
 */
template<int dim>
void TransportProblem<dim>::output_results()
{
   // output grid
   //------------
   if (parameters.output_meshes) {
      if (dim > 1)
         output_grid();
   }

   // create output data object
   //--------------------------
   DataOut<dim> data_out;
   data_out.attach_dof_handler(dof_handler);
   data_out.add_data_vector(present_solution, "flux");
   data_out.build_patches(degree + 1);

   // create output filename
   //-----------------------
   std::string output_extension;
   if (dim == 1) output_extension = ".gpl";
   else          output_extension = ".vtk";
   std::string viscosity_string;
   switch (parameters.viscosity_type) {
      case 0: {
         viscosity_string = "none";
         break;
      } case 1: {
         viscosity_string = "first_order";
         break;
      } case 2: {
         viscosity_string = "entropy";
         break;
      } case 3: {
         viscosity_string = "max_principle";
         break;
      } default: {
         Assert(false, ExcNotImplemented());
         break;
      }
   }
   std::stringstream output_filename_ss;
   output_filename_ss << "output/solution_" << viscosity_string;
   if (parameters.exact_solution_id != 0) output_filename_ss << "_" << parameters.exact_solution_id;
   output_filename_ss << output_extension;
   std::string output_filename = output_filename_ss.str();
   char *output_filename_char = (char*)output_filename.c_str();
   std::ofstream output(output_filename_char);

   // write solution output file
   //---------------------------
   if (dim == 1) data_out.write_gnuplot(output);
   else          data_out.write_vtk(output);

   // write viscosity output file
   //----------------------------
   if (parameters.viscosity_type != 0) {
      DataOut<dim> visc_out;
      visc_out.attach_dof_handler(dof_handler);
      // add viscosity data vector(s)
      if ((parameters.viscosity_type == 1)||(parameters.viscosity_type == 2))
         visc_out.add_data_vector(max_viscosity,"Max_Viscosity",DataOut<dim>::type_cell_data);
      if (parameters.viscosity_type == 2)
         visc_out.add_data_vector(entropy_viscosity,"Entropy_Viscosity",DataOut<dim>::type_cell_data);
      if (parameters.viscosity_type == 3)
         visc_out.add_data_vector(max_principle_viscosity,"Max_Principle_Viscosity",DataOut<dim>::type_cell_data);
      // build patches and write to file
      visc_out.build_patches(degree + 1);
      if (dim == 1) {
         std::ofstream visc_out_stream("output/viscosity.gpl");
         visc_out.write_gnuplot(visc_out_stream);
      } else {
         std::ofstream visc_out_stream("output/viscosity.vtk");
         visc_out.write_vtk(visc_out_stream);
      }
   }

   // print convergence table
   //------------------------
   if (parameters.exact_solution_id != 0) {
      convergence_table.set_precision("cell size", 3);
      convergence_table.set_scientific("cell size", true);
      convergence_table.set_precision("L2 error", 3);
      convergence_table.set_scientific("L2 error", true);
      convergence_table.evaluate_convergence_rates("L2 error", ConvergenceTable::reduction_rate_log2);
      std::cout << std::endl;
      convergence_table.write_text(std::cout);
   }
}

/** \brief evaluate error between numerical and exact solution
 */
template<int dim>
void TransportProblem<dim>::evaluate_error(const unsigned int cycle)
{
   // error per cell
   Vector<double> difference_per_cell (triangulation.n_active_cells());

   // compute error with analytic solution
   switch (parameters.exact_solution_id) {
      case 1: { // Test problem 1
         ExactSolution1<dim> exact_solution;

         VectorTools::integrate_difference (MappingQ<dim>(1),
	                                    dof_handler,
	                                    present_solution,
	                                    exact_solution,
	                                    difference_per_cell,
	                                    QGauss<dim>(degree+1),
	                                    VectorTools::L2_norm);
         break;
      } case 2: { // Test problem 2
         ExactSolution2<dim> exact_solution;

         VectorTools::integrate_difference (MappingQ<dim>(1),
                                            dof_handler,
                                            present_solution,
                                            exact_solution,
                                            difference_per_cell,
                                            QGauss<dim>(degree+1),
                                            VectorTools::L2_norm);
         break;
      } case 4: { // Test problem 4
         ExactSolution4<dim> exact_solution;

         VectorTools::integrate_difference (MappingQ<dim>(1),
                                            dof_handler,
                                            present_solution,
                                            exact_solution,
                                            difference_per_cell,
                                            QGauss<dim>(degree+1),
                                            VectorTools::L2_norm);
         break;
      } default: {
         Assert(false, ExcNotImplemented())
         break;
      }
   }

   // compute L2 error of vector of cell errors
   const double L2_error = difference_per_cell.l2_norm();

   const unsigned int n_active_cells = triangulation.n_active_cells();
   const unsigned int n_dofs = dof_handler.n_dofs();
   const double avg_cell_length = std::pow(2.0,dim) / std::pow(n_active_cells,1.0/dim);

   // add error values to convergence table
   convergence_table.add_value("cycle", cycle);
   convergence_table.add_value("cells", n_active_cells);
   convergence_table.add_value("dofs", n_dofs);
   convergence_table.add_value("cell size", avg_cell_length);
   convergence_table.add_value("L2 error", L2_error);
}

/** \brief check that the solution is non-negative at all nodes.
 */
template<int dim>
void TransportProblem<dim>::check_solution_nonnegative() const
{
   // loop over all degrees of freedom
   bool solution_is_negative = false;
   const unsigned int n_dofs = dof_handler.n_dofs();
   for (unsigned int i = 0; i < n_dofs; ++i)
      if (present_solution(i) < 0)
         solution_is_negative = true;

   // report if solution was negative anywhere or not
   if (solution_is_negative)
      std::cout << "Solution is negative!" << std::endl;
   else
      std::cout << "Solution is not negative at any node." << std::endl;
}

/** \brief check that the local discrete max principle is satisfied.
 */
template<int dim>
void TransportProblem<dim>::check_local_discrete_max_principle() const
{
   const unsigned int n_cells = triangulation.n_active_cells();
   Vector<double> max_values(n_cells, -1.0e15); // max of neighbors for each dof
   Vector<double> min_values(n_cells,  1.0e15); // min of neighbors for each dof

   std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
   for (; cell != endc; ++cell) {
      // get local dof indices
      cell->get_dof_indices(local_dof_indices);

      // find min and max values on cell
      double max_cell = -1.0e15; // initialized to arbitrary small value
      double min_cell =  1.0e15; // initialized to arbitrary large value
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

   // check that each dof value is bounded by its neighbors
   bool local_max_principle_satisfied = true;
   for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) {
      double value_i = present_solution(local_dof_indices[i]);
      if (value_i < min_values(i))
         local_max_principle_satisfied = false;
      if (value_i > max_values(i))
         local_max_principle_satisfied = false;
   }
   
   // report if local discrete maximum principle was satisfied or not
   if (local_max_principle_satisfied)
      std::cout << "The local discrete maximum principle was satisfied at all degrees of freedom." << std::endl;
   else
      std::cout << "The local discrete maximum principle was NOT satisfied at all degrees of freedom." << std::endl;
}

