#ifndef RefinementHandler_cc
#define RefinementHandler_cc

#include <deal.II/grid/tria.h>

using namespace dealii;

/** \brief Class for performing spatial/temporal refinement.
 */
template<int dim>
class RefinementHandler {
   public:

      enum RefinementMode { space, time };

      RefinementHandler(
         Triangulation<dim>   &triangulation,
         unsigned int         &n_cells,
         double               &dt,
         const RefinementMode refinement_mode,
         const bool           use_adaptive_refinement,
         const double         time_refinement_factor);
      ~RefinementHandler();

      void refine(unsigned int cycle) const;

   private:

      void refine_grid() const;

      Triangulation<dim> * const triangulation;
      unsigned int * const       n_cells;
      double * const             dt;
      const RefinementMode       refinement_mode;
      const bool                 use_adaptive_refinement;
      const double               time_refinement_factor;
};

#include "RefinementHandler.cc"
#endif
