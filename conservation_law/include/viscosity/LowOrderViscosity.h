/**
 * \file LowOrderViscosity.h
 * \brief Provides the header for the LowOrderViscosity class.
 */
#ifndef LowOrderViscosity_h
#define LowOrderViscosity_h

#include "include/viscosity/Viscosity.h"
#include "include/viscosity/ViscosityMultiplier.h"

using namespace dealii;

/**
 * \brief Class for low-order viscosity.
 *
 * The low-order viscosity is computed as
 * \f[
 *   \nu_K^L = c_{max} h_K \lambda_{K,max} \,,
 * \f]
 * where \f$\lambda_{K,max}\f$ is the maximum flux speed on cell \f$K\f$.
 * This viscosity is to be used in a standard Laplacian diffusion term:
 * \f[
 *   u_t + \nabla\cdot \mathbf{f}(u) - \nabla(\nu^L\nabla u) = 0 \,,
 * \f]
 * which makes the following contribution to the steady-state residual for
 * degree of freedom \f$i\f$, \f$r_i\f$:
 * \f[
 *   r_i = r_i + \int\limits_{S_i}\varphi_i\nabla(\nu\nabla u)dV \,.
 * \f]
 */
template <int dim>
class LowOrderViscosity : public Viscosity<dim>
{
public:
  /** \brief Alias for cell iterator */
  using Cell = typename Viscosity<dim>::Cell;

  /** \brief Alias for cell iterator map to double */
  using CellMap = typename Viscosity<dim>::CellMap;

  LowOrderViscosity(
    const double & c_max,
    CellMap & max_flux_speed,
    const FESystem<dim> & fe,
    const DoFHandler<dim> & dof_handler,
    const QGauss<dim> & cell_quadrature,
    const std::shared_ptr<ViscosityMultiplier<dim>> & viscosity_multiplier);

  void update(const Vector<double> & new_solution,
              const Vector<double> & old_solution,
              const double & dt,
              const unsigned int & n) override;

private:
  const double c_max;

  CellMap * const max_flux_speed;

  /** \brief Pointer to finite element system */
  const FESystem<dim> * fe;

  /** \brief Pointer to cell quadrature */
  const QGauss<dim> * cell_quadrature;

  /** \brief Pointer to viscosity multiplier */
  std::shared_ptr<ViscosityMultiplier<dim>> viscosity_multiplier;
};

#include "src/viscosity/LowOrderViscosity.cc"

#endif
