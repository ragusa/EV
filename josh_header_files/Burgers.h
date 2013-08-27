/** \file Burgers.h
 *  \brief Provides the header for the Burgers class.
 */
#ifndef Burgers_h
#define Burgers_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_component_interpretation.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/lac/vector.h>

#include "ConservationLaw.h"
#include "BurgersParameters.h"

/** \class Burgers
 *  \brief Provides the necessary functions for Burgers' equation.
 *
 *  This class extends the ConservationLaw class to viscous
 *  Burgers' equation:
 *  \f[
 *    u_t + u u_x = \nu u_{xx}
 *  \f]
 */
template <int dim>
class Burgers : public ConservationLaw<dim>
{
  public:
    Burgers(const BurgersParameters<dim> &params);

  private:
    BurgersParameters<dim> burgers_parameters;

    // number of components and position of components in solution vector
    static const unsigned int n_burgers_components = 1;

    const FEValuesExtractors::Scalar velocity;

    // vector of names of each component
    std::vector<std::string> get_component_names();

    // data component interpretation (scalar or vector) for outputting solution
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
       get_component_interpretations();

    void define_problem();

    void compute_cell_ss_residual(FEValues<dim> &fe_values,
                                  const typename DoFHandler<dim>::active_cell_iterator &cell,
                                  Vector<double> &cell_residual);
    void compute_face_ss_residual(FEFaceValues<dim> &fe_face_values,
                                  const typename DoFHandler<dim>::active_cell_iterator &cell,
                                  Vector<double> &cell_residual);
    Tensor<1,dim> flux_derivative(const double u);
    double entropy           (const double u) const;
    double entropy_derivative(const double u) const;
};

#include "Burgers.cc"

#endif
