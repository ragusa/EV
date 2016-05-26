/**
 * \file BoundaryDistance.h
 * \brief Provides the header for the BoundaryDistance class.
 */
#ifndef BoundaryDistance_cc
#define BoundaryDistance_cc

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <include/other/Exceptions.h>

using namespace dealii;

/**
 * \brief Base class for computing distance to the boundary along a line.
 */
template <int dim>
class BoundaryDistance
{
public:
  BoundaryDistance(const Tensor<1, dim> & direction);

  virtual double compute(const Point<dim> & x) const = 0;

protected:
  /** \brief transport direction */
  const Tensor<1, dim> direction;
};

#include "src/geometry/BoundaryDistance.cc"
#endif
