#ifndef RefinementHandler_cc
#define RefinementHandler_cc

#include <deal.II/grid/tria.h>
#include "TransportParameters.h"

using namespace dealii;

/**
 * Class for performing spatial/temporal refinement.
 */
template <int dim>
class RefinementHandler
{
public:
  RefinementHandler(const TransportParameters<dim> & parameters,
                    Triangulation<dim> & triangulation);

  void refine(unsigned int cycle);

  double get_nominal_time_step_size() const;

private:
  void refineGrid() const;

  /** pointer to triangulation */
  Triangulation<dim> * const triangulation;

  /** \brief Option to refine space in each refinement cycle */
  const bool refine_space;

  /** \brief Option to refine time in each refinement cycle */
  const bool refine_time;

  /** flag to use adaptive mesh refinement */
  const bool use_adaptive_refinement;

  /** factor to be applied to time step size if it is refined */
  const double time_refinement_factor;

  /** nominal time step size for current cycle */
  double nominal_dt;
};

#include "RefinementHandler.cc"
#endif
