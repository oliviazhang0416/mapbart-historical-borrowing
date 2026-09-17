#ifndef LRCBART_DRMU_H
#define LRCBART_DRMU_H

#include <algorithm>

// ---------------------------------------------------------------------------
// LRC-BART MODIFICATION START
// This keeps mBART's bottom-node enumeration and observation-to-leaf map, but
// replaces heterdrmu()'s two source means and per-leaf tau draw with the final
// scalar f draw or blocked (z,theta) g draw.
// ---------------------------------------------------------------------------
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi,
             double sigma1_sq, double sigma2_sq,
             int& K0, int& Lg, double& S0,
             rn& gen)
{
  tree::npv bnv;
  t.getbots(bnv);
  const size_t nb=bnv.size();

  std::vector<size_t> n1(nb,0), n2(nb,0);
  std::vector<double> sum1(nb,0.0), sum2(nb,0.0);
  std::map<tree::tree_p,size_t> bnmap;
  for(size_t l=0;l<nb;l++) bnmap[bnv[l]]=l;

  for(size_t i=0;i<di.n;i++) {
    double *xx=di.x+i*di.p;
    const size_t l=bnmap[t.bn(xx,xi)];
    if(di.source[i]==1) {
      ++n1[l];
      sum1[l]+=di.y[i];
    } else if(di.source[i]==2 && pi.kind==F_ENSEMBLE) {
      ++n2[l];
      sum2[l]+=di.y[i];
    }
  }

  for(size_t l=0;l<nb;l++) {
    if(pi.kind==F_ENSEMBLE) {
      const double A=(double)n1[l]/sigma1_sq+(double)n2[l]/sigma2_sq;
      const double B=sum1[l]/sigma1_sq+sum2[l]/sigma2_sq;
      const double denom=1.0+pi.lambda_f_sq*A;
      const double mean=pi.lambda_f_sq*B/denom;
      const double variance=pi.lambda_f_sq/denom;
      bnv[l]->settheta(mean+std::sqrt(variance)*gen.normal());
      bnv[l]->setz(1);
      continue;
    }

    // g leaf: draw z first, then theta conditional on z.
    const double pspike=g_pspike(
      n1[l],sum1[l],sigma1_sq,pi.tau0_sq,pi.tau1_sq,pi.w);
    const int z=(gen.uniform()<pspike) ? 1 : 0;
    const double tau_sq=z ? pi.tau0_sq : pi.tau1_sq;
    double theta;
    if(n1[l]==0) {
      // Single-arm/prior-only leaf: its posterior is exactly its prior.
      theta=std::sqrt(tau_sq)*gen.normal();
    } else {
      const double ybar=sum1[l]/(double)n1[l];
      const double v1=sigma1_sq/(double)n1[l];
      const double mean=tau_sq*ybar/(tau_sq+v1);
      const double variance=tau_sq*v1/(tau_sq+v1);
      theta=mean+std::sqrt(variance)*gen.normal();
    }
    bnv[l]->settheta(theta);
    bnv[l]->setz(z);

    ++Lg;
    if(z==1) { ++K0; S0+=theta*theta; }
  }
}
// LRC-BART MODIFICATION END

#endif
