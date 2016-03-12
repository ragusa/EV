/**
 * \file DMPExplicitEulerFCTFilter.h
 * \brief Provides the header for the DMPExplicitEulerFCTFilter class.
 */
#ifndef DMPExplicitEulerFCTFilter_h
#define DMPExplicitEulerFCTFilter_h

using namespace dealii;

/**
 * \brief Class for discrete maximum principle (DMP) filter for explicit Euler.
 */
template <int dim>
class DMPExplicitEulerFCTFilter : ExplicitEulerFCTFilter
{
public:
  DMPExplicitEulerFCTFilter();
};

#include "src/fct/DMPExplicitEulerFCTFilter.cc"

#endif
