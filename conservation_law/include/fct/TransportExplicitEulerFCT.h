/**
 * \file TransportExplicitEulerFCT.h
 * \brief Provides the header for the TransportExplicitEulerFCT class.
 */
#ifndef TransportExplicitEulerFCT_h
#define TransportExplicitEulerFCT_h

#include "include/parameters/TransportRunParameters.h"

using namespace dealii;

/**
 * \brief Class for ...
 */
template <int dim>
class TransportExplicitEulerFCT
{
public:
  TransportExplicitEulerFCT();

protected:
  virtual create_filters(const std::string & filter_sequence_string);
};

#include "src/fct/TransportExplicitEulerFCT.cc"

#endif
