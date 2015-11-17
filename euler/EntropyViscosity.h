/**
 * \file EntropyViscosity.h
 * \brief Provides the header for the EntropyViscosity class.
 */
#ifndef EntropyViscosity_h
#define EntropyViscosity_h

#include "Entropy.h"
#include "Viscosity.h"

using namespace dealii;

/**
 * \class EntropyViscosity
 * \brief Class for entropy viscosity.
 */
template <int dim>
class EntropyViscosity : public Viscosity<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Viscosity<dim>::Cell;

  EntropyViscosity(const DoFHandler<dim> & dof_handler);

  void update() override;

private:
  std::shared_ptr<Entropy<dim>> entropy;
};

#include "EntropyViscosity.cc"

#endif
