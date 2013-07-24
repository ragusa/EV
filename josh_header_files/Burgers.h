#ifndef Burgers_h
#define Burgers_h

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/vector.h>

#include "ConservationLaw.h"
#include "BurgersParameters.h"

/** \class Burgers
 *  \brief Provides the necessary functions for Burgers equation
 */
template <int dim>
class Burgers : public ConservationLaw<dim>
{
  public:
    Burgers(ParameterHandler &prm);

  private:
    BurgersParameters<dim> burgers_parameters;

    // number of components and position of components in solution vector
    static const unsigned int n_burgers_components = 1;

/*
    // vector of names of each component
    std::vector<std::string>
    get_component_names ();

    // data component interpretation (scalar or vector component) for outputting solution
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_component_interpretations ();
*/

    // compute flux matrix f(c)
    void compute_flux_matrix (const Vector<double> &W,
                              double (&flux)[n_euler_components][dim]);

    // computes forcing vector functions g(c)
    void compute_forcing_vector (const Vector<double> &W,
                                 double (&forcing)[n_euler_components]);

};

#include "Burgers.cc"

#endif
