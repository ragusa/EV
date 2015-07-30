#ifndef RefinementHandler_cc
#define RefinementHandler_cc

#include <deal.II/grid/tria.h>
#include "TransportParameters.h"

using namespace dealii;

/**
 * Class for performing spatial/temporal refinement.
 */
template<int dim>
class RefinementHandler
{
public:

  RefinementHandler(
    const TransportParameters<dim> & parameters,
    Triangulation<dim> & triangulation);
  ~RefinementHandler();

  void refine(unsigned int cycle) const;
  //double getNominalTimeStepSize() const;

private:

  void refineGrid() const;

  Triangulation<dim> * const triangulation;
  //double dt_nominal;
  const typename TransportParameters<dim>::RefinementMode refinement_mode;
  const bool use_adaptive_refinement;
  const double time_refinement_factor;
};

#include "RefinementHandler.cc"
#endif
