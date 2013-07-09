#include <iostream>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>

using namespace dealii;

template <int dim>
ConservationLaw<dim>::ConservationLaw(unsigned int n_components):
   n_components(n_components),
   mapping(),
   fe(FE_Q<dim>(1), n_components),
   dof_handler(triangulation),
   quadrature(2),
   face_quadrature(2),
   verbose_cout(std::cout, false)
{}

template <int dim>
void ConservationLaw<dim>::run()
{
   // make grid and refine
   GridGenerator::hyper_cube(triangulation,-1,1);
   triangulation.refine_global(3);

   // clear and distribute dofs
   dof_handler.clear();
   dof_handler.distribute_dofs(fe);

   // resize vectors
   old_old_solution.reinit(dof_handler.n_dofs());
   old_solution.reinit(dof_handler.n_dofs());
   current_solution.reinit(dof_handler.n_dofs());
   right_hand_side.reinit(dof_handler.n_dofs());

   // setup system
   setup_system();

   assemble_system();
//   std::pair<unsigned int, double> solution = solve();
   output_results();
}

template <int dim>
void ConservationLaw<dim>::setup_system ()
{
/*
   CompressedSparsityPattern sparsity_pattern (dof_handler.n_dofs(),
                                               dof_handler.n_dofs());
   DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
 
   system_matrix.reinit (sparsity_pattern);
*/
   std::cout << "Setup system" << std::endl;
}

template <int dim>
void ConservationLaw<dim>::assemble_system ()
{
   std::cout << "Assemble system" << std::endl;
}

template <int dim>
void ConservationLaw<dim>::assemble_cell_term (const FEValues<dim>             &fe_v,
                                               const std::vector<unsigned int> &dofs)
{}

template <int dim>
void ConservationLaw<dim>::assemble_face_term (const unsigned int               face_no,
                                               const FEFaceValuesBase<dim>     &fe_v,
                                               const FEFaceValuesBase<dim>     &fe_v_neighbor,
                                               const std::vector<unsigned int> &dofs,
                                               const std::vector<unsigned int> &dofs_neighbor,
                                               const bool                       external_face,
                                               const unsigned int               boundary_id,
                                               const double                     face_diameter)
{}

template <int dim>
std::pair<unsigned int, double> ConservationLaw<dim>::solve (Vector<double> &solution)
{
   std::pair<unsigned int, double> solutions;
   std::cout << "Solve" << std::endl;
   return solutions;
}

template <int dim>
void ConservationLaw<dim>::compute_refinement_indicators (Vector<double> &indicator) const
{}

template <int dim>
void ConservationLaw<dim>::refine_grid (const Vector<double> &indicator)
{}

template <int dim>
void ConservationLaw<dim>::output_results () const
{
   std::cout << "Output results" << std::endl;
}
