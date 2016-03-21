#ifndef TransportThetaExecutioner_cc
#define TransportThetaExecutioner_cc

#include "include/executioners/TransportTransientExecutioner.h"
#include "include/fct/TransportThetaFCT.h"
#include "include/parameters/RunParameters.h"

using namespace dealii;

/**
 * \brief Class for theta executioner of a transport problem.
 */
template <int dim>
class TransportThetaExecutioner : public TransportTransientExecutioner<dim>
{
public:
  /** \brief Alias for scheme */
  using Scheme = typename TransportExecutioner<dim>::Scheme;
  /** \brief Alias for high-order scheme */
  using HighOrderScheme = typename TransportExecutioner<dim>::HighOrderScheme;

  TransportThetaExecutioner(const TransportRunParameters & parameters,
                            TransportProblemParameters<dim> & problem_parameters,
                            Triangulation<dim> & triangulation,
                            PostProcessor<dim> & postprocessor,
                            const double & nominal_dt);

protected:
  void compute_new_solution(const double & dt,
                            const double & dt_old,
                            const double & t_old,
                            const unsigned int & n) override;

  void compute_galerkin_solution(const double & dt, const double & t_new);

  void compute_low_order_solution(const double & dt, const double & t_new);

  void compute_high_order_solution(const double & dt,
                                   const double & t_new,
                                   const unsigned int & n);

  void compute_entropy_viscosity_solution(const double & dt,
                                          const double & t_new,
                                          const unsigned int & n);

  void compute_fct_solution(const double & dt, const double & t_old);

  /** \brief new high-order diffusion matrix \f$\mathbf{D}^{H,n+1}\f$ */
  SparseMatrix<double> high_order_diffusion_matrix_new;

  /** \brief new high-order steady-state matrix \f$\mathbf{A}^{H,n+1}\f$ */
  SparseMatrix<double> high_order_ss_matrix_new;

  /** \brief theta parameter \f$\theta\f$ */
  const double theta;

private:
  void perform_fct_ssprk_step(const double & dt,
                              const double & old_stage_dt,
                              const unsigned int & n);

  /** \brief FCT */
  std::shared_ptr<TransportThetaFCT<dim>> fct;
};

#include "src/executioners/TransportThetaExecutioner.cc"
#endif
