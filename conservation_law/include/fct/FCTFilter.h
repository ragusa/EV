/**
 * \file FCTFilter.h
 * \brief Provides the header for the FCTFilter class.
 */
#ifndef FCTFilter_h
#define FCTFilter_h

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include "include/fct/DoFBounds.h"
#include "include/fct/Limiter.h"
#include "include/parameters/RunParameters.h"

using namespace dealii;

/**
 * \brief Abstract base class for FCT filters, which are used to limit
 *        antidiffusion fluxes.
 */
template <int dim>
class FCTFilter
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  FCTFilter(const RunParameters & run_parameters,
            const std::shared_ptr<Limiter<dim>> limiter,
            const DoFHandler<dim> & dof_handler,
            const FESystem<dim> & fe);

  virtual bool check_bounds(const Vector<double> & new_solution) const;

  Vector<double> get_lower_solution_bound() const;

  Vector<double> get_upper_solution_bound() const;

protected:
  void compute_min_and_max_of_dof_vector(const Vector<double> & dof_vector,
                                         Vector<double> & min_values,
                                         Vector<double> & max_values) const;

  void enforce_antidiffusion_bounds_signs();

  /** \brief solution bounds \f$\mathbf{W}^\pm\f$ */
  DoFBounds<dim> solution_bounds;

  /** \brief antidiffusion bounds \f$\mathbf{Q}^\pm\f$ */
  DoFBounds<dim> antidiffusion_bounds;

  /** \brief limiter */
  const std::shared_ptr<Limiter<dim>> limiter;

  /** \brief degree of freedom handler */
  const DoFHandler<dim> * const dof_handler;

  /** \brief number of degrees of freedom */
  const unsigned int n_dofs;

  /** \brief number of components */
  const unsigned int n_components;

  /** \brief number of degrees of freedom per cell */
  const unsigned int dofs_per_cell;

  /** \brief number of degrees of freedom per cell per component */
  const unsigned int dofs_per_cell_per_component;

  /** \brief option to force correct signs of antidiffusion bounds */
  bool do_enforce_antidiffusion_bounds_signs;
};

#include "src/fct/FCTFilter.cc"

#endif
