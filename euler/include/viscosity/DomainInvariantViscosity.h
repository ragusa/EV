/**
 * \file DomainInvariantViscosity.h
 * \brief Provides the header for the DomainInvariantViscosity class.
 */
#ifndef DomainInvariantViscosity_h
#define DomainInvariantViscosity_h

using namespace dealii;

/**
 * \class DomainInvariantViscosity
 * \brief Class for ...
 */
template <int dim>
class DomainInvariantViscosity
{
public:
  DomainInvariantViscosity(const Triangulation<dim> & triangulation);

protected:
  double compute_max_wave_speed(const Tensor<1,dim> normal,
    const std::vector<double> solution_i, const std::vector<double> solution_j);

  /** \brief Degree of freedom handler for a scalar */
  DoFHandler<dim> dof_handler;

  /** \brief Pointer to triangulation */
  const Triangulation<dim> * const triangulation;

  /** \brief Number of lines in triangulation */
  unsigned int n_lines;

  /** \brief Vector of gradient tensors from 1st node on each line */
  std::vector<Tensor<1,dim>> c_ij;

  /** \brief Vector of gradient tensors from 2nd node on each line */
  std::vector<Tensor<1,dim>> c_ji;

  /** \brief Vector of norm of gradient tensors from 1st node on each line */
  std::vector<double> c_ij_norm;

  /** \brief Vector of norm of gradient tensors from 2nd node on each line */
  std::vector<double> c_ji_norm;

  /** \brief Vector of normal vectors from 1st node on each line */
  std::vector<Tensor<1,dim>> normal_ij;

  /** \brief Vector of normal vectors from 2nd node on each line */
  std::vector<Tensor<1,dim>> normal_ji;
};

#include "src/viscosity/DomainInvariantViscosity.cc"

#endif
