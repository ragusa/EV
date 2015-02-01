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
   const FunctionParser<dim> &source_function,
   const std::string     &entropy_string,
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
      entropy_residual_coefficient(entropy_residual_coefficient),
      jump_coefficient(jump_coefficient),
      domain_volume(domain_volume),
      entropy_viscosity(n_cells)
{
   // initialize entropy function
   std::map<std::string,double> constants;
   entropy_function.initialize("u",entropy_string,constants,false);

   // initialize transport direction
   //transport_direction = transport_direction;
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
void EntropyViscosity<dim>::compute_entropy_domain_average(const Vector<double> &new_solution,
                                                           const Vector<double> &old_solution)
{
   FEValues<dim> fe_values(*fe, cell_quadrature, update_values | update_JxW_values);

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
   max_entropy_deviation_domain = 0.0;

   // loop over cells
   for (cell = (*dof_handler).begin_active(); cell != endc; ++cell) {
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

/** \brief Computes the entropy viscosity for each cell.
 */
template <int dim>
Vector<double> EntropyViscosity<dim>::compute_entropy_viscosity(const Vector<double> &new_solution,
                                                                const Vector<double> &old_solution,
                                                                const double         &dt)
{
   // compute entropy average in domain
   compute_entropy_domain_average(new_solution,old_solution);

   std::vector<Point<dim> > points(n_q_points_cell);
   std::vector<double> cross_section_values(n_q_points_cell);
   std::vector<double> source_values(n_q_points_cell);

   // FE cell values for computing entropy
   FEValues<dim> fe_values(*fe, cell_quadrature,
         update_values | update_gradients | update_quadrature_points
               | update_JxW_values);

   // FE face values for computing entropy jumps
   FEFaceValues<dim> fe_values_face         (*fe, face_quadrature,
      update_values | update_gradients | update_JxW_values | update_normal_vectors);
   FEFaceValues<dim> fe_values_face_neighbor(*fe, face_quadrature,
      update_values | update_gradients | update_JxW_values | update_normal_vectors);

   // cell values
   std::vector<double>         new_values   (n_q_points_cell);
   std::vector<Tensor<1,dim> > new_gradients(n_q_points_cell);
   std::vector<double>         old_values   (n_q_points_cell);
   std::vector<double>         new_entropy_values(n_q_points_cell,0.0);
   std::vector<double>         old_entropy_values(n_q_points_cell,0.0);

   // face values
   std::vector<double>         values_face            (n_q_points_face);
   std::vector<double>         values_face_neighbor   (n_q_points_face);
   std::vector<Tensor<1,dim> > gradients_face         (n_q_points_face);
   std::vector<Tensor<1,dim> > gradients_face_neighbor(n_q_points_face);
   std::vector<Point<dim> >    normal_vectors         (n_q_points_face);
   Vector<double>              entropy_face           (n_q_points_face);

   // cell iterator
   typename DoFHandler<dim>::active_cell_iterator cell = (*dof_handler).begin_active(),
                                                  endc = (*dof_handler).end();
   // loop over cells
   unsigned int i_cell = 0;
   for (cell = (*dof_handler).begin_active(); cell != endc; ++cell, ++i_cell)
   {
      // reinitialize FE values
      fe_values.reinit(cell);

      // get previous time step (n) values and gradients
      fe_values[flux].get_function_values   (new_solution, new_values);
      fe_values[flux].get_function_gradients(new_solution, new_gradients);
      // get previous previous time step (n-1) values
      fe_values[flux].get_function_values   (old_solution, old_values);
   
      // get quadrature points on cell
      points = fe_values.get_quadrature_points();

      // get source values for all quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         source_values[q] = (*source_function).value(points[q]);

      // get cross section values for all quadrature points
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         cross_section_values[q] = (*cross_section_function).value(points[q]);

      // compute max entropy residual in cell
      //----------------------------------------------------------------------------
      // compute entropy values at each quadrature point on cell. The entropy definition s = 0.5*u^2 is used.
      for (unsigned int q = 0; q < n_q_points_cell; ++q) {
         new_entropy_values[q] = 0.5 * new_values[q] * new_values[q];
         old_entropy_values[q] = 0.5 * old_values[q] * old_values[q];
      }
      // compute entropy residual values at each quadrature point on cell
      std::vector<double> entropy_residual_values(n_q_points_cell,0.0);
      for (unsigned int q = 0; q < n_q_points_cell; ++q)
         entropy_residual_values[q] = (new_entropy_values[q] - old_entropy_values[q])/dt
            + transport_direction * new_values[q] * new_gradients[q]
            + cross_section_values[q] * new_values[q] * new_values[q]
            - source_values[q] * new_values[q];
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
      entropy_viscosity(i_cell) = (entropy_residual_coefficient * max_entropy_residual
         + jump_coefficient * max_jump_in_cell) / max_entropy_deviation_domain;
   }

   return entropy_viscosity;
}
