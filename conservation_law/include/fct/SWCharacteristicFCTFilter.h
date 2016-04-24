/**
 * \file SWCharacteristicFCTFilter.h
 * \brief Provides the header for the SWCharacteristicFCTFilter class.
 */
#ifndef SWCharacteristicFCTFilter_h
#define SWCharacteristicFCTFilter_h

#include "include/fct/CharacteristicFCTFilter.h"

using namespace dealii;

/**
 * \brief Shallow water characteristic FCT filter.
 */
template <int dim>
class SWCharacteristicFCTFilter : public CharacteristicFCTFilter<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  SWCharacteristicFCTFilter(
    const RunParameters & run_parameters,
    const std::shared_ptr<Limiter<dim>> limiter,
    const DoFHandler<dim> & dof_handler,
    const FESystem<dim> & fe,
    const SparseMatrix<double> & lumped_mass_matrix,
    const double & gravity,
    const std::map<unsigned int, double> & dirichlet_values);

protected:
  FullMatrix<double> compute_transformation_matrix(
    const Vector<double> & solution) const override;

  /** \brief Acceleration due to gravity \f$g\f$ */
  const double gravity;
};

#include "src/fct/SWCharacteristicFCTFilter.cc"

#endif
