/**
 * \file HyperboxBoundaryDistance.h
 * \brief Provides the header for the HyperboxBoundaryDistance class.
 */
#ifndef HyperboxBoundaryDistance_cc
#define HyperboxBoundaryDistance_cc

#include "include/geometry/BoundaryDistance.h"

using namespace dealii;

/**
 * \brief Base class for computing distance to the boundary along a line.
 */
template <int dim>
class HyperboxBoundaryDistance : public BoundaryDistance<dim>
{
public:
  HyperboxBoundaryDistance(const Tensor<1, dim> & direction,
                           const double & x_min,
                           const double & y_min,
                           const double & z_min);

  double compute(const Point<dim> & x) const override;

protected:
  /** \brief lower boundaries of the hyperbox in each dimension */
  std::vector<double> x_min;
};

#include "src/geometry/HyperboxBoundaryDistance.cc"
#endif
