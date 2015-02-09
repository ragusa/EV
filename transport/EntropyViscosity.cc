/** \brief constructor
 */
template<int dim>
EntropyViscosity<dim>::EntropyViscosity(
   const FESystem<dim>   &fe,
   const unsigned int    &n_cells,
   const DoFHandler<dim> &dof_handler,
   const QGauss<dim>     &cell_quadrature,
   const QGauss<dim-1>   &face_quadrature,
   const Tensor<1,dim>   &transport_direction,
   const FunctionParser<dim> &cross_section_function,
   FunctionParser<dim>   &source_function,
   const std::string     &entropy_string,
   const std::string     &entropy_derivative_string,
   const double          &entropy_residual_coefficient,
   const double          &jump_coefficient,
   const double          &domain_volume) :                      
      fe(&fe),
      flux(0),
      n_cells(n_cells),
      dof_handler(&dof_handler),
      n_dofs(dof_handler.n_dofs()),
      dofs_per_cell(fe.dofs_per_cell),
      faces_per_cell(GeometryInfo<dim>::faces_per_cell),
      cell_quadrature(cell_quadrature),
      face_quadrature(face_quadrature),
      n_q_points_cell(cell_quadrature.size()),
      n_q_points_face(face_quadrature.size()),
      transport_direction(transport_direction),
      cross_section_function(&cross_section_function),
      source_function(&source_function),
      entropy_string(entropy_string),
      entropy_derivative_string(entropy_derivative_string),
      entropy_residual_coefficient(entropy_residual_coefficient),
      jump_coefficient(jump_coefficient),
      domain_volume(domain_volume),
      entropy_viscosity(n_cells)
{
   // initialize entropy function
   std::map<std::string,double> constants;
   entropy_function.initialize("u",entropy_string,constants,false);
   entropy_derivative_function.initialize("u",entropy_derivative_string,constants,false);
}

/** \brief destructor
 */
template<int dim>
EntropyViscosity<dim>::~EntropyViscosity() {
}

/** \brief computes the domain-averaged entropy and the max entropy
 *         deviation in the domain
 */
template <int dim>
void EntropyViscosity<dim>::compute_entropy_domain_average(const Vector<double> &old_solution)
{
   FEValues<dim> fe_values(*fe, cell_quadrature, update_values | update_JxW_values);
   std::vector<double>    old_solution_values      (n_q_points_cell);
   std::vector<double>    old_entropy              (n_q_points_cell);
   std::vector<Point<1> > old_solution_local_points(n_q_points_cell);

   // compute domain-averaged entropy
   //--------------------------------
   double domain_integral_entropy = 0.0;

   // loop over cells
   typename DoFHandler<dim>::active_cell_iterator cell = (*dof_handler).begin_active(),
                                                  endc = (*dof_handler).end();
   for (cell = (*dof_handler).begin_active(); cell != endc; ++cell) {
      // reinitialize FE values
      fe_values.reinit(cell);
      // get solution values
      fe_values[flux].get_function_values(old_solution,old_solution_values);
      // compute entropy values
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         old_solution_local_points[q][0] = old_solution_values[q];
      entropy_function.value_list(old_solution_local_points, old_entropy);
      // loop over quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         // add contribution of quadrature point to entropy integral
         domain_integral_entropy += old_entropy[q] * fe_values.JxW(q);
      }
   }
   // domain-averaged entropy
   domain_averaged_entropy = domain_integral_entropy / domain_volume;

   // compute max deviation of entropy from domain-averaged entropy
   //--------------------------------------------------------------
   max_entropy_deviation_domain = 0.0;

   // loop over cells
   for (cell = (*dof_handler).begin_active(); cell != endc; ++cell) {
      // reinitialize FE values
      fe_values.reinit(cell);
      // get old values
      std::vector<double> old_solution_values(n_q_points_cell);
      fe_values[flux].get_function_values(old_solution,old_solution_values);
      // compute entropy values
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         old_solution_local_points[q][0] = old_solution_values[q];
      entropy_function.value_list(old_solution_local_points, old_entropy);
      // loop over quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         // add contribution of quadrature point to entropy integral
         max_entropy_deviation_domain = std::max(max_entropy_deviation_domain,
               std::abs(old_entropy[q] - domain_averaged_entropy));
      }
   }
}

/** \brief Computes the entropy viscosity for each cell.
 */
template <int dim>
Vector<double> EntropyViscosity<dim>::compute_entropy_viscosity(const Vector<double> &old_solution,
                                                                const Vector<double> &older_solution,
                                                                const double         &dt,
                                                                const double         &time)
{
   // set the time of the source function
   source_function->set_time(time);
   
   // compute entropy average in domain
   compute_entropy_domain_average(old_solution);

   // cell values
   std::vector<Point<dim> > points              (n_q_points_cell);
   std::vector<double>      cross_section_values(n_q_points_cell);
   std::vector<double>      source_values       (n_q_points_cell);
   std::vector<double>      old_solution_local  (n_q_points_cell);
   std::vector<double>      older_solution_local(n_q_points_cell);
   std::vector<double>      old_entropy         (n_q_points_cell);
   std::vector<double>      older_entropy       (n_q_points_cell);
   std::vector<double>      old_entropy_derivative(n_q_points_cell);
   std::vector<Point<1> >   old_solution_local_points  (n_q_points_cell);
   std::vector<Point<1> >   older_solution_local_points(n_q_points_cell);
   std::vector<Tensor<1,dim> > old_gradients(n_q_points_cell);

   // face values
   std::vector<double>         values_face            (n_q_points_face);
   std::vector<double>         values_face_neighbor   (n_q_points_face);
   std::vector<Tensor<1,dim> > gradients_face         (n_q_points_face);
   std::vector<Tensor<1,dim> > gradients_face_neighbor(n_q_points_face);
   std::vector<Point<dim> >    normal_vectors         (n_q_points_face);
   Vector<double>              entropy_face           (n_q_points_face);

   // FE cell values for computing entropy
   FEValues<dim> fe_values(*fe, cell_quadrature,
         update_values | update_gradients | update_quadrature_points
               | update_JxW_values);

   // FE face values for computing entropy jumps
   FEFaceValues<dim> fe_values_face         (*fe, face_quadrature,
      update_values | update_gradients | update_JxW_values | update_normal_vectors);
   FEFaceValues<dim> fe_values_face_neighbor(*fe, face_quadrature,
      update_values | update_gradients | update_JxW_values | update_normal_vectors);

   // cell iterator
   typename DoFHandler<dim>::active_cell_iterator cell = (*dof_handler).begin_active(),
                                                  endc = (*dof_handler).end();
   // loop over cells
   unsigned int i_cell = 0;
   for (cell = (*dof_handler).begin_active(); cell != endc; ++cell, ++i_cell)
   {
      // reinitialize FE values
      fe_values.reinit(cell);

      // get solution values and gradients
      fe_values[flux].get_function_values   (old_solution,   old_solution_local);
      fe_values[flux].get_function_values   (older_solution, older_solution_local);
      fe_values[flux].get_function_gradients(old_solution,   old_gradients);
   
      // get cross section and source values for all quadrature points
      points = fe_values.get_quadrature_points();
      source_function       ->value_list(points,source_values);
      cross_section_function->value_list(points,cross_section_values);

      // compute max entropy residual in cell
      //----------------------------------------------------------------------------
      // compute entropy values at each quadrature point on cell
      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         old_solution_local_points[q][0]   = old_solution_local[q];
         older_solution_local_points[q][0] = older_solution_local[q];
      }
      entropy_function.value_list(old_solution_local_points,  old_entropy);
      entropy_function.value_list(older_solution_local_points,older_entropy);
      entropy_derivative_function.value_list(old_solution_local_points,old_entropy_derivative);

      // compute entropy residual values at each quadrature point on cell
      std::vector<double> entropy_residual_values(n_q_points_cell,0.0);
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         entropy_residual_values[q] = (old_entropy[q] - older_entropy[q])/dt
            + old_entropy_derivative[q] * (transport_direction * old_gradients[q]
            + cross_section_values[q] * old_solution_local[q]
            - source_values[q]);

      // determine maximum entropy residual in cell
      double max_entropy_residual = 0.0;
      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         max_entropy_residual  = std::max(max_entropy_residual,  std::abs(entropy_residual_values[q]));
      }
   
      // compute max jump in cell
      //----------------------------------------------------------------------------
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
            fe_values_face.get_function_values(old_solution, values_face);
   
            // compute entropy at each quadrature point on face
            for (unsigned int q = 0; q < n_q_points_face; ++q)
               entropy_face[q] = 0.5 * values_face[q] * values_face[q];
   
            // get gradients on adjacent faces of current cell and neighboring cell
            fe_values_face.get_function_gradients(         old_solution, gradients_face);
            fe_values_face_neighbor.get_function_gradients(old_solution, gradients_face_neighbor);
   
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
      entropy_viscosity(i_cell) = (entropy_residual_coefficient * max_entropy_residual
         + jump_coefficient * max_jump_in_cell) / max_entropy_deviation_domain;
   }

   return entropy_viscosity;
}
