/**
 * \file FCTFilter.h
 * \brief Provides the header for the FCTFilter class.
 */
#ifndef FCTFilter_h
#define FCTFilter_h

#include "include/fct/DoFBounds.h"

using namespace dealii;

/**
 * \brief Abstract base class for FCT filters, which are used to limit
 *        antidiffusion fluxes.
 */
template <int dim>
class FCTFilter : public DoFBounds<dim>
{
public:
  FCTFilter();
};

#include "src/fct/FCTFilter.cc"

#endif
