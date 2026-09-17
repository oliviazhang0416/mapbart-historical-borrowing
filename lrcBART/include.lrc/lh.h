#ifndef LRCBART_LH_H
#define LRCBART_LH_H

#include <limits>

// ---------------------------------------------------------------------------
// LRC-BART MODIFICATION START
// This replaces mBART's heterlh() shared-leaf marginal with the two exact
// final-revision marginals.  All calculations are in log space.
// ---------------------------------------------------------------------------

double ldnorm(double x, double mean, double variance)
{
  static const double LOG_2PI=1.8378770664093453;
  return -0.5*(LOG_2PI+std::log(variance)+(x-mean)*(x-mean)/variance);
}

double lse2(double a, double b)
{
  const double neg_inf=-std::numeric_limits<double>::infinity();
  if(a==neg_inf) return b;
  if(b==neg_inf) return a;
  const double mx=std::max(a,b);
  return mx+std::log(std::exp(a-mx)+std::exp(b-mx));
}

// Heteroscedastic f-leaf marginal.  Constants common to all tree structures
// cancel in the Metropolis ratio, leaving the Chipman factor below.
double f_lh(size_t n1, size_t n2,
                     double s1, double s2,
                     double sigma1_sq, double sigma2_sq,
                     double lambda_f_sq)
{
  const double A=(double)n1/sigma1_sq+(double)n2/sigma2_sq;
  const double B=s1/sigma1_sq+s2/sigma2_sq;
  const double denom=1.0+lambda_f_sq*A;
  return -0.5*std::log(denom)+lambda_f_sq*B*B/(2.0*denom);
}

// Spike/slab g-leaf marginal relative to the no-leaf-effect density.
double g_lh(size_t n1, double s1, double sigma1_sq,
                     double tau0_sq, double tau1_sq, double w)
{
  if(n1==0) return 0.0;
  const double ybar=s1/(double)n1;
  const double v1=sigma1_sq/(double)n1;
  const double neg_inf=-std::numeric_limits<double>::infinity();
  const double l0=(w>0.0) ? std::log(w)+ldnorm(ybar,0.0,tau0_sq+v1) : neg_inf;
  const double l1=(w<1.0) ? std::log(1.0-w)+ldnorm(ybar,0.0,tau1_sq+v1) : neg_inf;
  return lse2(l0,l1)-ldnorm(ybar,0.0,v1);
}

// Posterior probability that a g leaf is in the spike state z=1.
double g_pspike(size_t n1, double s1, double sigma1_sq,
                         double tau0_sq, double tau1_sq, double w)
{
  if(w>=1.0) return 1.0;
  if(w<=0.0) return 0.0;
  if(n1==0) return w;
  const double ybar=s1/(double)n1;
  const double v1=sigma1_sq/(double)n1;
  const double l0=std::log(w)+ldnorm(ybar,0.0,tau0_sq+v1);
  const double l1=std::log(1.0-w)+ldnorm(ybar,0.0,tau1_sq+v1);
  return std::exp(l0-lse2(l0,l1));
}

double lh(size_t n1, size_t n2,
                   double s1, double s2,
                   double sigma1_sq, double sigma2_sq,
                   const pinfo& pi)
{
  if(pi.kind==G_ENSEMBLE)
    return g_lh(n1,s1,sigma1_sq,pi.tau0_sq,pi.tau1_sq,pi.w);
  return f_lh(n1,n2,s1,s2,sigma1_sq,sigma2_sq,pi.lambda_f_sq);
}

// LRC-BART MODIFICATION END

#endif
