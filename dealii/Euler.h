/** \file Euler.h
 *  \brief Provides the header for the Euler class.
 */
#ifndef Euler_h
#define Euler_h

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
#include "EulerParameters.h"

/** \class Euler
 *  \brief Provides the necessary functions for Euler equations.
 */
template <int dim>
class Euler : public ConservationLaw<dim>
{
  public:
    Euler(const EulerParameters<dim> &params);

  private:
    EulerParameters<dim> euler_parameters;

    // number of components and position of components in solution vector
    static const unsigned int n_euler_components = dim+2;

    const FEValuesExtractors::Scalar density;
    const FEValuesExtractors::Vector velocity;
    const FEValuesExtractors::Scalar energy;

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

#include "Euler.cc"

#endif
