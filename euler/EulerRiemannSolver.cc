/**
 * \file EulerRiemannSolver.cc
 * \brief Provides the function definitions for the EulerRiemannSolver class.
 */

/**
 * \brief Constructor.
 */
template<int dim>
EulerRiemannSolver<dim>::EulerRiemannSolver(
  const double &rho_left_,
  const double &u_left_,
  const double &p_left_,
  const double &rho_right_,
  const double &u_right_,
  const double &p_right_,
  const double &gamma_,
  const double &x_interface_
) :
    Function<dim>(dim+2),
    rho_left(rho_left_),
    u_left(u_left_),
    p_left(p_left_),
    rho_right(rho_right_),
    u_right(u_right_),
    p_right(p_right_),
    gamma(gamma_),
    x_interface(x_interface_),
    g1((gamma-1)/(2*gamma)),
    g2((gamma+1)/(2*gamma)),
    g3(2*gamma/(gamma-1)),
    g4(2/(gamma-1)),
    g5(2/(gamma+1)),
    g6((gamma-1)/(gamma+1)),
    g7((gamma-1)/2),
    g8(gamma-1),
    c_left(std::sqrt(gamma*p_left/rho_left)),
    c_right(std::sqrt(gamma*p_right/rho_right))
{
  // check for vacuum conditions
  Assert(g4*(c_left+c_right)>(u_right-u_left), ExcInvalidState())

  // get exact value of the pressure in the star region
  // first, find educated guess of the value of the pressure in the star region
  const double cup = 0.25*(rho_left+rho_right)*(c_left+c_right);
  double ppv = 0.5*(p_left+p_right)+0.5*(u_left-u_right)*cup;
  ppv = std::max(0.0,ppv);
  const double pmin = std::min(p_left,p_right);
  const double pmax = std::max(p_left,p_right);
  const double qmax = pmax/pmin;
  const double quser = 2;
  if ((qmax <= quser) && ((pmin <= ppv) && (ppv <= pmax)))
  { // select PVRS Riemann solver
    p_star = ppv;
  }
  else
  {
    if (ppv < pmin)
    { // select two-rarefaction Riemann solver
      const double pq = std::pow((p_left/p_right),g1);
      u_star = (pq*u_left/c_left+u_right/c_right+g4*(pq-1))/
        (pq/c_left+1/c_right);
      const double ptl = 1+g7*(u_left-u_star)/c_left;
      const double ptr = 1+g7*(u_star-u_right)/c_right;
      p_star = 0.5*(p_left*std::pow(ptl,g3) + p_right*std::pow(ptr,g3));
    }
    else
    { // select two-shock Riemann solver with PVRS as estimate
      const double gel = std::sqrt((g5/rho_left)/(g6*p_left+ppv));
      const double ger = std::sqrt((g5/rho_right)/(g6*p_right+ppv));
      p_star = (gel*p_left+ger*p_right-(u_right-u_left))/(gel+ger);
    }
  }

  // second, perform Newton iteration
  const unsigned int max_iter = 5000;
  const double tolpre = 1.0e-12;
  double pold = p_star;
  bool converged;
  double fl;
  double fr;
  for (unsigned int i = 0; i < max_iter; ++i)
  {
    // evaluate the pressure function on the left and on the right
    double fld;
    if (pold <= p_left)
    { // rarefaction
        double pratl = pold/p_left;
        fl = g4*c_left*(std::pow(pratl,g1)-1);
        fld = (1/(rho_left*c_left))*std::pow(pratl,-g2);
    }
    else
    { // shock
        double al = g5/rho_left;
        double bl = g6*p_left;
        double qrtl = std::sqrt(al/(bl+pold));
        fl = (pold-p_left)*qrtl;
        fld = (1-0.5*(pold-p_left)/(bl+pold))*qrtl;
    }

    double frd;
    if (pold <= p_right)
    { // rarefaction
        double pratr = pold/p_right;
        fr = g4*c_right*(std::pow(pratr,g1)-1);
        frd = (1/(rho_right*c_right))*std::pow(pratr,-g2);
    }
    else
    { // shock
        double ar = g5/rho_right;
        double br = g6*p_right;
        double qrtr = std::sqrt(ar/(br+pold));
        fr = (pold-p_right)*qrtr;
        frd = (1-0.5*(pold-p_right)/(br+pold))*qrtr;
    }

    p_star = pold-(fl+fr+u_right-u_left)/(fld+frd);
    const double change = 2*std::abs((p_star-pold)/(p_star+pold));
    if (change <= tolpre)
    {
      converged = true;
      break;
    } else {
      converged = false;
    }
    if (p_star < 0)
      p_star = tolpre;

    pold = p_star;
  }
  
  // make sure that iteration converged
  Assert(converged, ExcInvalidState());
  
  // find the velocity in the star region
  u_star = 0.5*(u_left+u_right + fr-fl);
}

/** Computes exact solution at a single point.
 */
template<int dim>
double EulerRiemannSolver<dim>::value(
  const Point<dim>   &p,
  const unsigned int component) const
{
  // get time
  double t = this->get_time();

  // get x position
  double x = p[0];

  // sample the solution
  double rho_sample;
  double u_sample;
  double p_sample;
  if (t == 0)
  {
    if (x <= x_interface)
    {
      rho_sample = rho_left;
      u_sample = u_left;
      p_sample = p_left;
    }
    else
    {
      rho_sample = rho_right;
      u_sample = u_right;
      p_sample = p_right;
    }
  }
  else
  {
    double s = (x-x_interface)/t;
    if (s <= u_star) // left of contact discontinuity
    {
      if (p_star <= p_left) // left rarefaction
      {
        double shl = u_left-c_left;
        if (s <= shl) // left data state
        {
          rho_sample = rho_left;
          u_sample = u_left;
          p_sample = p_left;
        }
        else
        {
          double cml = c_left*std::pow(p_star/p_left,g1);
          double stl = u_star-cml;
          if (s > stl) // star left state
          {
            rho_sample = rho_left*std::pow(p_star/p_left,1/gamma);
            u_sample = u_star;
            p_sample = p_star;
          }
          else // left fan
          {
            u_sample = g5*(c_left+g7*u_left+s);
            double cs = g5*(c_left+g7*(u_left-s));
            rho_sample = rho_left*std::pow(cs/c_left,g4);
            p_sample = p_left*std::pow(cs/c_left,g3);
          }
        }
      }
      else // left shock
      {
        double pml = p_star/p_left;
        double sl = u_left-c_left*std::sqrt(g2*pml+g1);
        if (s <= sl) // left data state
        {
          rho_sample = rho_left;
          u_sample = u_left;
          p_sample = p_left;
        }
        else // left star state
        {
          rho_sample = rho_left*(pml+g6)/(pml*g6+1);
          u_sample = u_star;
          p_sample = p_star;
        }
      }
    }
    else // right of contact discontinuity
    {
      if (p_star > p_right) // right shock
      {
        double pmr = p_star/p_right;
        double sr = u_right+c_right*std::sqrt(g2*pmr+g1);
        if (s >= sr) // right data state
        {
          rho_sample = rho_right;
          u_sample = u_right;
          p_sample = p_right;
        }
        else // star right state
        {
          rho_sample = rho_right*(pmr+g6)/(pmr*g6+1);
          u_sample = u_star;
          p_sample = p_star;
        }
      }
      else // right rarefaction
      {
        double shr = u_right+c_right;
        if (s >= shr) // right data state
        {
          rho_sample = rho_right;
          u_sample = u_right;
          p_sample = p_right;
        }
        else
        {
          double cmr = c_right*std::pow(p_star/p_right,g1);
          double str = u_star+cmr;
          if (s <= str) // star right state
          {
            rho_sample = rho_right*std::pow(p_star/p_right,1/gamma);
            u_sample = u_star;
            p_sample = p_star;
          }
          else // right fan
          {
            u_sample = g5*(-c_right+g7*u_right+s);
            double cs = g5*(c_right-g7*(u_right-s));
            rho_sample = rho_right*std::pow(cs/c_right,g4);
            p_sample = p_right*std::pow(cs/c_right,g3);
          }
        }
      }
    }
  }

  // return solution for the requested component
  if (component == 0)
    return rho_sample;
  else if (component == 1)
    return rho_sample*u_sample;
  else
  {
    const double e_sample = p_sample/((gamma-1.0)*rho_sample);
    const double E_sample = (e_sample + 0.5*std::pow(u_sample,2))*rho_sample;
    return E_sample;
  }
}
