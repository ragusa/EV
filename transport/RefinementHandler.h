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
  ~RefinementHandler();

  void refine(unsigned int cycle) const;

private:
  void refineGrid() const;

  /** pointer to triangulation */
  Triangulation<dim> * const triangulation;

  /** type of refinement to occur in each cycle, either spatial or temporal */
  const typename TransportParameters<dim>::RefinementMode refinement_mode;

  /** flag to use adaptive mesh refinement */
  const bool use_adaptive_refinement;

  /** factor to be applied to time step size if it is refined */
  const double time_refinement_factor;
};

#include "RefinementHandler.cc"
#endif
