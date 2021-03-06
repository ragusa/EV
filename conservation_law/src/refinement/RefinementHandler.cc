/**
 * \brief Constructor.
 */
template <int dim>
// RefinementHandler<dim>::RefinementHandler(const RunParameters & parameters_,
RefinementHandler<dim>::RefinementHandler(const RunParameters & parameters_,
                                          Triangulation<dim> & triangulation_)
  : triangulation(&triangulation_),
    refine_space(parameters_.refine_space),
    refine_time(parameters_.refine_time),
    use_adaptive_refinement(parameters_.use_adaptive_mesh_refinement),
    time_refinement_factor(parameters_.time_refinement_factor),
    nominal_dt(parameters_.time_step_size)
{
}

/**
 * \brief Refines the grid or time step size.
 *
 * \param[in] cycle  cycle index, first is 0
 */
template <int dim>
void RefinementHandler<dim>::refine(unsigned int cycle)
{
  if (cycle != 0)
  {
    // refine mesh
    if (refine_space)
      refineGrid();

    // refine time
    if (refine_time)
      nominal_dt *= time_refinement_factor;
  }
}

/**
 * \brief Refines the grid.
 */
template <int dim>
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

/**
 * \brief Returns the nominal time step size for this cycle
 *
 * \return nominal time step size for this cycle
 */
template <int dim>
double RefinementHandler<dim>::get_nominal_time_step_size() const
{
  return nominal_dt;
}
