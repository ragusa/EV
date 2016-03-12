/**
 * \file DoFBounds.h
 * \brief Provides the header for the DoFBounds class.
 */
#ifndef DoFBounds_h
#define DoFBounds_h

using namespace dealii;

/**
 * \brief Abstract base class for implementing upper and lower bounds for
 *        some degree of freedom indexed vector such as a solution vector.
 */
class DoFBounds
{
public:
  DoFBounds(const unsigned int & n_dofs);

  double get_min(const unsigned int & i) const;

  double get_max(const unsigned int & i) const;

  void widen(const DoFBounds & dof_bounds);

  void shrink(const DoFBounds & dof_bounds);

  void remove_bounds(const std::vector<double> & dof_indices);

protected:
  /** \brief lower bound vector */
  Vector<double> lower_bound;

  /** \brief upper bound vector */
  Vector<double> upper_bound;

  /** \brief number of degrees of freedom */
  const unsigned int n_dofs;
};

#include "src/fct/DoFBounds.cc"

#endif
