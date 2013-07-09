#include <iostream>

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
   std::cout << "Running..." << std::endl;
}

template <int dim>
void ConservationLaw<dim>::setup_system ()
{}

template <int dim>
void ConservationLaw<dim>::assemble_system ()
{}

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
{}
