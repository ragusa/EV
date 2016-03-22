/**
 * \file ShallowWaterExplicitEulerFCT.h
 * \brief Provides the header for the ShallowWaterExplicitEulerFCT class.
 */
#ifndef ShallowWaterExplicitEulerFCT_h
#define ShallowWaterExplicitEulerFCT_h

#include "include/fct/ExplicitEulerFCT.h"
#include "include/fct/SWCharacteristicFCTFilter.h"
#include "include/parameters/ShallowWaterRunParameters.h"
#include "include/parameters/ShallowWaterProblemParameters.h"

using namespace dealii;

/**
 * \brief Class for explicit Euler FCT for the shallow water equations
 */
template <int dim>
class ShallowWaterExplicitEulerFCT : public ExplicitEulerFCT<dim>
{
public:
  ShallowWaterExplicitEulerFCT(
    const ShallowWaterRunParameters & run_parameters,
    const ShallowWaterProblemParameters<dim> & problem_parameters,
    const DoFHandler<dim> & dof_handler_,
    const FESystem<dim> & fe,
    const SparseMatrix<double> & consistent_mass_matrix,
    const SparseMatrix<double> & lumped_mass_matrix);

protected:
  virtual std::shared_ptr<ExplicitEulerFCTFilter<dim>> create_filter(
    const std::string & filter_string) override;

  /** \brief Acceleration due to gravity \f$g\f$ */
  const double gravity;
};

#include "src/fct/ShallowWaterExplicitEulerFCT.cc"

#endif
