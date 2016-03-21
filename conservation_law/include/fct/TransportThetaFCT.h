/**
 * \file TransportThetaFCT.h
 * \brief Provides the header for the TransportThetaFCT class.
 */
#ifndef TransportThetaFCT_h
#define TransportThetaFCT_h

#include "include/fct/ThetaFCT.h"
#include "include/parameters/TransportRunParameters.h"
#include "include/parameters/TransportProblemParameters.h"

using namespace dealii;

/**
 * \brief Class for theta FCT for linear transport.
 */
template <int dim>
class TransportThetaFCT : public ThetaFCT<dim>
{
public:
  TransportThetaFCT(const TransportRunParameters & run_parameters,
                    const TransportProblemParameters<dim> & problem_parameters,
                    const DoFHandler<dim> & dof_handler_,
                    const FESystem<dim> & fe,
                    const SparseMatrix<double> & consistent_mass_matrix,
                    const SparseMatrix<double> & lumped_mass_matrix);

protected:
  virtual std::shared_ptr<ThetaFCTFilter<dim>> create_filter(
    const std::string & filter_string) override;
};

#include "src/fct/TransportThetaFCT.cc"

#endif
