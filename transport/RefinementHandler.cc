/**
 * Constructor.
 */
template<int dim>
RefinementHandler<dim>::RefinementHandler(
  const TransportParameters<dim> & parameters_,
  Triangulation<dim> & triangulation_) :
    triangulation(& triangulation_),
    refinement_mode(parameters_.refinement_mode),
    use_adaptive_refinement(parameters_.use_adaptive_refinement),
    time_refinement_factor(parameters_.time_refinement_factor)
{
}

/**
 * Destructor.
 */
template<int dim>
RefinementHandler<dim>::~RefinementHandler()
{
}

/**
 * Refine the grid or time step size.
 */
template<int dim>
void RefinementHandler<dim>::refine(unsigned int cycle) const
{
  if (cycle != 0)
  {
    // refine mesh if user selected the option
    switch (refinement_mode)
    {
      case TransportParameters<dim>::RefinementMode::space :
      { // refine space
        refineGrid();
        break;
      }
      case TransportParameters<dim>::RefinementMode::time :
      { // refine time
        ExcNotImplemented();
        //dt_nominal *= time_refinement_factor;
        break;
      }
      default :
      {
        ExcNotImplemented();
      }
    }
  }
}

/**
 * Refines the grid.
 */
template<int dim>
void RefinementHandler<dim>::refineGrid() const
{

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
}
