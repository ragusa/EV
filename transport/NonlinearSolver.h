#ifndef NonlinearSolver_cc
#define NonlinearSolver_cc

#include "LinearSolver.h"
#include "TransportParameters.h"

using namespace dealii;

/**
 * Class for nonlinear solver.
 */
template<int dim>
class NonlinearSolver
{
public:

  NonlinearSolver(const TransportParameters<dim> & parameters);

private:

  void endIteration();
  bool checkConvergence();

  const double nonlinear_tolerance;
  const double iteration_max;

  Vector<double> previous_solution;
  unsigned int iteration_number;
};

#include "NonlinearSolver.cc"
#endif
