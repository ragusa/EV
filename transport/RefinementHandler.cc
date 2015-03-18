/** \brief Constructor.
 */
template<int dim>
RefinementHandler<dim>::RefinementHandler(
   Triangulation<dim>   &triangulation,
   unsigned int         &n_cells,
   double               &dt,
   const RefinementMode refinement_mode,
   const bool           use_adaptive_refinement,
   const double         time_refinement_factor) :
   triangulation(&triangulation),
   n_cells(&n_cells),
   dt(&dt),
   refinement_mode(refinement_mode),
   use_adaptive_refinement(use_adaptive_refinement),
   time_refinement_factor(time_refinement_factor)
{
}

/** \brief Destructor.
 */
template<int dim>
RefinementHandler<dim>::~RefinementHandler()
{
}

/** \brief Refine the grid or time step size.
 */
template<int dim>
void RefinementHandler<dim>::refine(unsigned int cycle) const {
   if (cycle != 0) {
      // refine mesh if user selected the option
      switch (refinement_mode) {
         case RefinementMode::space : { // refine space
            refine_grid();
            break;
         } case RefinementMode::time : { // refine time
            *dt = *dt * time_refinement_factor;
            break;
         } default : {
            ExcNotImplemented();
         }
      }
   }
}

/** \brief Refine the grid.
 */
template<int dim>
void RefinementHandler<dim>::refine_grid() const {

   Assert(use_adaptive_refinement == false, ExcNotImplemented());
/*
   // refine adaptively or uniformly
   if (use_adaptive_mesh_refinement) { // adaptively
      Vector<float> estimated_error_per_cell(*n_cells);

      KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim-1>(3),
            typename FunctionMap<dim>::type(), new_solution,
            estimated_error_per_cell);

      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
            estimated_error_per_cell, 0.3, 0.03);

      triangulation.execute_coarsening_and_refinement();
   } else // uniformly
*/
   // refine uniformly
   triangulation->refine_global(1);

   // update number of cells
   *n_cells = triangulation->n_active_cells();
}
