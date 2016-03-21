/**
 * \file TransportExplicitEulerFCT.h
 * \brief Provides the header for the TransportExplicitEulerFCT class.
 */
#ifndef TransportExplicitEulerFCT_h
#define TransportExplicitEulerFCT_h

#include "include/fct/ExplicitEulerFCT.h"
#include "include/parameters/TransportRunParameters.h"
#include "include/parameters/TransportProblemParameters.h"

using namespace dealii;

/**
 * \brief Class for explicit Euler FCT for linear transport.
 */
template <int dim>
class TransportExplicitEulerFCT : public ExplicitEulerFCT<dim>
{
public:
  TransportExplicitEulerFCT(
    const TransportRunParameters & run_parameters,
    const TransportProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler_,
    const FESystem<dim> & fe,
    const SparseMatrix<double> & consistent_mass_matrix,
    const SparseMatrix<double> & lumped_mass_matrix);

protected:
  virtual std::shared_ptr<ExplicitEulerFCTFilter<dim>> create_filter(
    const std::string & filter_string) override;
};

#include "src/fct/TransportExplicitEulerFCT.cc"

#endif
