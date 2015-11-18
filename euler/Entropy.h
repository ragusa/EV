/**
 * \file Entropy.h
 * \brief Provides the header for the Entropy class.
 */
#ifndef Entropy_h
#define Entropy_h

using namespace dealii;

/**
 * \class Entropy
 * \brief Class for entropy.
 */
template <int dim>
class Entropy
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename DoFHandler<dim>::active_cell_iterator;

  Entropy();

  void reinitialize(const Vector<double> & solution);

  virtual std::vector<double> compute_entropy() const = 0;

  virtual std::vector<double> compute_divergence_entropy_flux() const = 0;

  virtual std::vector<double> compute_entropy_normalization() const;

private:
  virtual void reinitialize_group_fe_values(const Vector<double> & solution) {}
  const bool need_to_compute_average_entropy;
};

#include "Entropy.cc"

#endif
