/** \file Burgers.cc
 *  \brief Provides the function definitions for the Burgers class.
 */

/** \fn Burgers<dim>::Burgers(ParameterHandler &prm)
 *  \brief Constructor for the Burgers class.
 *
 *  In addition to initializing member data, this function gets
 *  Burgers parameters.
 */
template <int dim>
Burgers<dim>::Burgers(ParameterHandler &prm):
   ConservationLaw<dim>(prm,n_burgers_components)
{
   // get Burgers parameters
   burgers_parameters.get_parameters(prm);
   // initialize initial conditions
   this->initial_conditions.initialize (FunctionParser<dim>::default_variable_names(),
                                        burgers_parameters.initial_conditions_expressions,
                                        std::map<std::string, double>());
   this->component_names = this->get_component_names();
   this->component_interpretations = this->get_component_interpretations();
}

/** \fn std::vector<std::string> Burgers<dim>::get_component_names()
 *  \brief Returns the names of each component.
 *  \return vector of names of each component
 */
template <int dim>
std::vector<std::string> Burgers<dim>::get_component_names ()
{
   std::vector<std::string> names (1, "velocity");

   return names;
}

/** \fn std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *      Burgers<dim>::get_component_interpretations()
 *  \brief Returns the interpretations for each component.
 *
 *  This function returns the interpretation of each component,
 *  i.e., whether each component is a scalar or a component of
 *  a vector.
 *  \return data component interpretations
 */
template <int dim>
std::vector<DataComponentInterpretation::DataComponentInterpretation>
   Burgers<dim>::get_component_interpretations ()
{
   std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (1, DataComponentInterpretation::component_is_scalar);

   return data_component_interpretation;
} 

/** \fn Burgers<dim>::compute_ss_residual(double t, Vector<double> &solution)
 *  \brief Computes the steady-state residual for the Burgers' equation.
 *
 *  This function computes the steady-state residual, where component \f$i\f$ is
 *  \f[
 *    ({\psi}_i,\mathbf{g}(\mathbf{u})-\mathbf{f}(\mathbf{u})).
 *  \f]
 *  For the inviscid Burgers' equation, this is the following:
 *  \f[
 *    ({\psi}_i,-u u_x)
 *  \f]
 *  \param t time at which the steady-state residual is to be evaluated
 *  \param solution at which to evaluate the steady-state residual
 */
template <int dim>
void Burgers<dim>::compute_ss_residual(double t, Vector<double> &solution)
{
   const FEValuesExtractors::Scalar velocity (0);

   QGauss<dim> quadrature_formula(5);
   FEValues<dim> fe_values (this->fe, quadrature_formula,
                            update_values | update_gradients | update_JxW_values);

   const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
   const unsigned int n_q_points    = quadrature_formula.size();

   Vector<double> cell_residual(dofs_per_cell);

   std::vector<unsigned int> local_dof_indices (dofs_per_cell);

   // unit vector in x-direction; used for recovering du/dx from gradient
   Tensor<1,dim> unit_x_vector; unit_x_vector[0] = 1;

   typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler.begin_active(),
                                                  endc = this->dof_handler.end();
   // loop over cells
   for (; cell!=endc; ++cell)
   {
      cell_residual = 0;
      fe_values.reinit(cell);
      // get global DoF indices
      cell->get_dof_indices(local_dof_indices);

      // loop over test functions
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
         // loop over components in value term
         for (unsigned int j = 0; j < dofs_per_cell; ++j)
            // loop over components in gradient term
            for (unsigned int m = 0; m < dofs_per_cell; ++m)
               // loop over quadrature points
               for (unsigned int q = 0; q < n_q_points; ++q)
                  cell_residual(i) -= fe_values[velocity].value(i,q)
                                      *solution(local_dof_indices[j])*fe_values[velocity].value(j,q)
                                      *solution(local_dof_indices[m])*fe_values[velocity].gradient(m,q)
                                      *unit_x_vector
                                      *fe_values.JxW(q);

      // add artificial viscosity if requested
      if (burgers_parameters.viscosity_type != BurgersParameters<dim>::none)
      {
         // compute viscosity
         Vector<double> viscosity(n_q_points);
         switch (burgers_parameters.viscosity_type)
         {
            case BurgersParameters<dim>::constant:
            {
               for (unsigned int q = 0; q < n_q_points; ++q)
                  viscosity(q) = burgers_parameters.constant_viscosity_value;
               break;
            }
            default:
            {
               Assert(false,ExcNotImplemented());
               break;
            }
         }
         // loop over test functions
         for (unsigned int i = 0; i < dofs_per_cell; ++i)
            // loop over linear combination functions
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
               // loop over quadrature points
               for (unsigned int q = 0; q < n_q_points; ++q)
                  cell_residual(i) -= fe_values[velocity].gradient(i,q)
                                      *viscosity(q)
                                      *solution(local_dof_indices[j])*fe_values[velocity].gradient(j,q)
                                      *fe_values.JxW(q);
      }

      // aggregate local residual into global residual
      // loop over test functions
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
         this->ss_residual(local_dof_indices[i]) += cell_residual(i);
   }

   // apply boundary conditions: zero Dirichlet
   std::map<unsigned int, double> boundary_values;
   VectorTools::interpolate_boundary_values(this->dof_handler,
                                            0,
                                            ZeroFunction<dim>(),
                                            boundary_values);
   for (std::map<unsigned int, double>::iterator bv = boundary_values.begin(); bv != boundary_values.end(); ++bv)
      this->ss_residual(bv->first) = (bv->second);
}
