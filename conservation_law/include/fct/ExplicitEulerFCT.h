/**
 * \file ExplicitEulerFCT.h
 * \brief Provides the header for the ExplicitEulerFCT class.
 */
#ifndef ExplicitEulerFCT_h
#define ExplicitEulerFCT_h

using namespace dealii;

/**
 * \brief Class for ...
 */
template <int dim>
class ExplicitEulerFCT
{
public:
  ExplicitEulerFCT();

  void compute_antidiffusion_vector(
    const Vector<double> & high_order_solution,
    const Vector<double> & old_solution,
    const SparseMatrix<double> & lumped_mass_matrix,
    const SparseMatrix<double> & consistent_mass_matrix,
    const double & dt,
    const SparseMatrix<double> & low_order_diffusion_matrix,
    const SparseMatrix<double> & high_order_diffusion_matrix,
    Vector<double> & antidiffusion);
};

#include "src/fct/ExplicitEulerFCT.cc"

#endif
