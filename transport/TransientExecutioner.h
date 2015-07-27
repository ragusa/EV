#ifndef TransientExecutioner_cc
#define TransientExecutioner_cc

#include "Executioner.h"

using namespace dealii;

/**
 * Class for transient executioner.
 */
template<int dim>
class TransientExecutioner : public Executioner<dim>
{
public:
  TransientExecutioner(const TransportParameters<dim> & parameters,
    const Triangulation<dim> & triangulation,
    const Tensor<1, dim> & transport_direction,
    const FunctionParser<dim> & cross_section_function,
    FunctionParser<dim> & source_function, Function<dim> & incoming_function);
  ~TransientExecutioner();

  void run() override;
};

#include "TransientExecutioner.cc"
#endif
