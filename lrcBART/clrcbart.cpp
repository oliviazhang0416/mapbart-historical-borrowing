#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <ctime>
#include <cmath>
#include <algorithm>

// Create aliases retained from mBART/cmbart.cpp.
#define printf Rprintf
#define cout Rcpp::Rcout
typedef std::vector<double> v1d;
typedef std::vector<v1d> v2d;
typedef std::vector<v2d> v3d;

using std::endl;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Header order follows mBART/cmbart.cpp.  The directory name is the only
// structural rename; all LRC changes live in the corresponding mBART files.
#include "include.lrc/rn.h"
#include "include.lrc/tree.h"
#include "include.lrc/treefuns.h"
#include "include.lrc/info.h"
#include "include.lrc/getsuff.h"
#include "include.lrc/lh.h"
#include "include.lrc/drmu.h"
#include "include.lrc/bd.h"
#include "include.lrc/bart.h"
#include "include.lrc/rtnorm.h"

// ---------------------------------------------------------------------------
// LRC-BART MODIFICATION START
// Direct sourceCpp descendant of cmbart().  x and x_test retain mBART's p x n
// convention (observations are contiguous columns).  The revision's many new
// sampler settings are grouped in hyper/control lists rather than introducing
// a package interface.
// ---------------------------------------------------------------------------
// JOINT LRC-BART MODIFICATION START
// Trailing optional treatment input preserves original positional calls.
// [[Rcpp::export]]
Rcpp::List clrcbart(
    Rcpp::NumericMatrix x,          // p x n training covariates
    Rcpp::NumericVector y,          // centered outcome or centered log-time
    Rcpp::IntegerVector source,     // 1 = trial (controls only in legacy mode), 2 = external
    Rcpp::IntegerVector status,     // 1 = observed, 0 = censored
    Rcpp::NumericVector censor,     // centered log censoring bound
    Rcpp::NumericMatrix x_test,     // p x n_test; may have zero columns
    Rcpp::List cutpoints,
    Rcpp::List hyper,
    Rcpp::List control,
    Rcpp::Nullable<Rcpp::NumericVector> treatment = R_NilValue)
// JOINT LRC-BART MODIFICATION END
{
  const size_t p=(size_t)x.nrow();
  const size_t n=(size_t)x.ncol();
  const size_t np=(size_t)x_test.ncol();
  if((size_t)y.size()!=n || (size_t)source.size()!=n)
    Rcpp::stop("x, y, and source have incompatible sizes");
  if(np>0 && (size_t)x_test.nrow()!=p)
    Rcpp::stop("x_test must have the same number of predictor rows as x");
  if(p==0 || n==0) Rcpp::stop("x must contain predictors and observations");

  Rcpp::IntegerVector event;
  if(status.size()==0) event=Rcpp::IntegerVector(n,1);
  else {
    if((size_t)status.size()!=n) Rcpp::stop("status has the wrong length");
    event=Rcpp::clone(status);
  }
  Rcpp::NumericVector lower;
  if(censor.size()==0) lower=Rcpp::NumericVector(n,NA_REAL);
  else {
    if((size_t)censor.size()!=n) Rcpp::stop("censor has the wrong length");
    lower=Rcpp::clone(censor);
  }

  size_t n1=0,n2=0;
  for(size_t i=0;i<n;i++) {
    if(source[i]==1) ++n1;
    else if(source[i]==2) ++n2;
    else Rcpp::stop("source values must be 1 (RCT) or 2 (RWD)");
  }
  // JOINT LRC-BART MODIFICATION START
  // PR#2 R/lrcbart.R: source codes describe membership, not treatment.
  // Only absence of ALL trial rows invokes the inherited prior-only g path.
  const bool no_trial_data=(n1==0);
  const bool joint_model=treatment.isNotNull();
  Rcpp::NumericVector A;
  size_t n_treated=0;
  double tau=0.0; // manuscript psi; distinct from discrepancy tau0/tau1
  double tau_prior_var=100.0;
  if(joint_model) {
    A=Rcpp::NumericVector(treatment.get());
    if((size_t)A.size()!=n)
      Rcpp::stop("treatment must have length ncol(x)");
    if(hyper.containsElementNamed("tau_prior_var"))
      tau_prior_var=Rcpp::as<double>(hyper["tau_prior_var"]);
    if(!std::isfinite(tau_prior_var) || tau_prior_var<=0.0)
      Rcpp::stop("tau_prior_var must be positive and finite");
    for(size_t i=0;i<n;i++) {
      if(!std::isfinite(A[i]) || (A[i]!=0.0 && A[i]!=1.0))
        Rcpp::stop("treatment must contain only finite 0/1 indicators");
      if(source[i]==2 && A[i]!=0.0)
        Rcpp::stop("external observations must have treatment zero");
      if(A[i]==1.0) ++n_treated;
      if(!std::isfinite(y[i])) Rcpp::stop("joint outcomes must be finite");
    }
  }
  // Joint single-arm means treated trial data without concurrent controls.
  // Preserve the old flag's meaning for callers omitting treatment.
  const bool single_arm=joint_model ? (n_treated>0 && n_treated==n1) : no_trial_data;
  // JOINT LRC-BART MODIFICATION END

  const size_t Hf=(size_t)Rcpp::as<int>(hyper["H_f"]);
  const size_t Hg=(size_t)Rcpp::as<int>(hyper["H_g"]);
  const double alpha_f=Rcpp::as<double>(hyper["alpha_f"]);
  const double beta_f=Rcpp::as<double>(hyper["beta_f"]);
  const double alpha_g=Rcpp::as<double>(hyper["alpha_g"]);
  const double beta_g=Rcpp::as<double>(hyper["beta_g"]);
  const double lambda_f_sq=Rcpp::as<double>(hyper["lambda_f_sq"]);
  const double nu_sigma=Rcpp::as<double>(hyper["nu_sigma"]);
  const double lambda_sigma1=Rcpp::as<double>(hyper["lambda_sigma1"]);
  const double lambda_sigma2=Rcpp::as<double>(hyper["lambda_sigma2"]);
  const double nu0=Rcpp::as<double>(hyper["nu0"]);
  const double s0_sq=Rcpp::as<double>(hyper["s0_sq"]);
  const double tau1_sq=Rcpp::as<double>(hyper["tau1_sq"]);
  const double aw=Rcpp::as<double>(hyper["a_w"]);
  const double bw=Rcpp::as<double>(hyper["b_w"]);
  const bool update_tau0=Rcpp::as<bool>(hyper["update_tau0"]);
  const bool augmentation=Rcpp::as<bool>(hyper["augmentation"]);
  const double p_grow=Rcpp::as<double>(hyper["p_grow"]);
  const double p_prune=Rcpp::as<double>(hyper["p_prune"]);
  const size_t n_min_f=(size_t)Rcpp::as<int>(hyper["n_min"]);
  const size_t n_min_g_requested=hyper.containsElementNamed("n_min_g") ?
    (size_t)Rcpp::as<int>(hyper["n_min_g"]) : n_min_f;
  // JOINT LRC-BART MODIFICATION START
  const size_t n_min_g=no_trial_data ? 0 : n_min_g_requested;
  // JOINT LRC-BART MODIFICATION END

  if(Hf==0) Rcpp::stop("H_f must be positive");
  if(lambda_f_sq<=0.0 || tau1_sq<=0.0 || s0_sq<=0.0)
    Rcpp::stop("leaf-prior variances must be positive");
  if(p_grow<0.0 || p_prune<0.0 || p_grow+p_prune>1.0)
    Rcpp::stop("invalid grow/prune probabilities");

  double sigma1_sq=Rcpp::as<double>(hyper["sigma1_sq_init"]);
  double sigma2_sq=Rcpp::as<double>(hyper["sigma2_sq_init"]);
  double tau0_sq=std::min(s0_sq,0.5*tau1_sq);
  const double wfix=Rcpp::as<double>(hyper["w_fixed"]);
  const bool w_fixed=!std::isnan(wfix);
  if(w_fixed && !std::isfinite(wfix))
    Rcpp::stop("w_fixed must be NA or finite");
  double w=w_fixed ? wfix : aw/(aw+bw);
  if(w<0.0 || w>1.0) Rcpp::stop("w_fixed must be NA or lie in [0,1]");

  const int n_burn=Rcpp::as<int>(control["n_burn"]);
  const int n_draw=Rcpp::as<int>(control["n_draw"]);
  const int thin=Rcpp::as<int>(control["thin"]);
  const bool verbose=Rcpp::as<bool>(control["verbose"]);
  const bool keep_trees=Rcpp::as<bool>(control["keep_trees"]);
  const bool keep_train=Rcpp::as<bool>(control["keep_train"]);
  const bool keep_test=Rcpp::as<bool>(control["keep_test"]);
  const int seed=Rcpp::as<int>(control["seed"]);
  if(n_burn<0 || n_draw<1 || thin<1) Rcpp::stop("invalid MCMC control values");

  // JOINT LRC-BART ADDITION START
  // PR#2 sampler control: repeat g updates, then update hyperparameters once.
  const double gs=control.containsElementNamed("g_sweeps") ?
    Rcpp::as<double>(control["g_sweeps"]) : 1.0;
  if(!std::isfinite(gs) || gs<1.0 || gs!=std::floor(gs) ||
     gs>(double)std::numeric_limits<int>::max())
    Rcpp::stop("g_sweeps must be a positive integer");
  const int g_sweeps=(int)gs;
  if(joint_model) {
    if(!std::isfinite(sigma1_sq) || sigma1_sq<=0.0 ||
       !std::isfinite(sigma2_sq) || sigma2_sq<=0.0)
      Rcpp::stop("joint initial source variances must be positive and finite");
    for(size_t i=0;i<n;i++) {
      if(event[i]!=0 && event[i]!=1)
        Rcpp::stop("status values must be 0 or 1");
      if(!augmentation && event[i]==0)
        Rcpp::stop("censored joint outcomes require augmentation=TRUE");
    }
  }
  // JOINT LRC-BART ADDITION END

  Rcpp::NumericVector latent=Rcpp::clone(y);
  if(augmentation) {
    for(size_t i=0;i<n;i++) {
      if(event[i]!=0 && event[i]!=1) Rcpp::stop("status values must be 0 or 1");
      if(event[i]==0) {
        if(!std::isfinite(lower[i])) Rcpp::stop("censored rows need finite censor bounds");
        const double sd=std::sqrt(source[i]==1 ? sigma1_sq : sigma2_sq);
        latent[i]=lower[i]+0.1*sd;
      }
    }
  }

  // Inherited cmbart.cpp pattern: copy the R cutpoint object directly into
  // the existing xinfo representation at sampler setup.
  if((size_t)cutpoints.size()!=p)
    Rcpp::stop("cutpoints must contain one numeric vector per predictor");
  xinfo xi(p);
  for(size_t v=0;v<p;v++) {
    Rcpp::NumericVector cv=cutpoints[v];
    xi[v].assign(cv.begin(),cv.end());
    if(xi[v].empty()) Rcpp::stop("each predictor needs at least one cutpoint");
  }
  std::vector<int> nc(p);
  for(size_t v=0;v<p;v++) nc[v]=(int)xi[v].size();

  bart fforest(Hf);
  bart gforest(Hg);
  fforest.setxinfo(xi);
  gforest.setxinfo(xi);
  fforest.setprior(F_ENSEMBLE,alpha_f,beta_f,lambda_f_sq,
                   tau0_sq,tau1_sq,w,n_min_f,p_grow,p_prune);
  gforest.setprior(G_ENSEMBLE,alpha_g,beta_g,lambda_f_sq,
                   tau0_sq,tau1_sq,w,n_min_g,p_grow,p_prune);
  fforest.setdata(p,n,&x[0],&latent[0],&source[0],&nc[0]);
  gforest.setdata(p,n,&x[0],&latent[0],&source[0],&nc[0]);

  Rcpp::NumericVector out_sig1(n_draw),out_sig2(n_draw),out_tau0(n_draw),out_w(n_draw);
  Rcpp::IntegerVector out_K0(n_draw),out_Lg(n_draw);
  // JOINT LRC-BART ADDITION START
  Rcpp::NumericVector out_tau(joint_model ? n_draw : 0);
  // JOINT LRC-BART ADDITION END
  Rcpp::NumericMatrix f_train(keep_train ? n_draw : 0,keep_train ? n : 0);
  Rcpp::NumericMatrix g_train(keep_train ? n_draw : 0,keep_train ? n : 0);
  Rcpp::NumericMatrix f_test(keep_test ? n_draw : 0,keep_test ? np : 0);
  Rcpp::NumericMatrix g_test(keep_test ? n_draw : 0,keep_test ? np : 0);
  Rcpp::List f_draws(keep_trees ? n_draw : 0);
  Rcpp::List g_draws((keep_trees && Hg>0) ? n_draw : 0);

  std::vector<double> offset_f(n,0.0),offset_g(n,0.0);
  std::vector<double> pred_f(np,0.0),pred_g(np,0.0);
  arn gen;
  gen.set_seed((unsigned int)seed);

  const int total=n_burn+n_draw*thin;
  int saved=0;
  for(int it=0;it<total;it++) {
    // JOINT LRC-BART MODIFICATION START
    // PR#2 compute_resid: all groups enter f after source/treatment adjustment.
    for(size_t i=0;i<n;i++) {
      offset_f[i]=(source[i]==1) ? gforest.f(i) : 0.0;
      if(joint_model) offset_f[i]+=A[i]*tau;
    }
    // JOINT LRC-BART MODIFICATION END
    int dummyK=0,dummyL=0;
    double dummyS0=0.0;
    fforest.draw(sigma1_sq,sigma2_sq,&offset_f[0],
                 dummyK,dummyL,dummyS0,gen);

    int K0=0,Lg=0;
    double S0=0.0;
    if(Hg>0) {
      // JOINT LRC-BART MODIFICATION START
      // Both trial arms enter g; the forest excludes source-2 rows as before.
      for(size_t i=0;i<n;i++) {
        offset_g[i]=fforest.f(i);
        if(joint_model) offset_g[i]+=A[i]*tau;
      }
      gforest.set_g_hyper(tau0_sq,tau1_sq,w);
      for(int sweep=0;sweep<g_sweeps;sweep++) {
        // drmu accumulates per-leaf statistics. Only the FINAL sweep's
        // retained leaves may contribute to the once-per-iteration update.
        K0=0; Lg=0; S0=0.0;
        gforest.draw(sigma1_sq,sigma2_sq,&offset_g[0],K0,Lg,S0,gen);
      }
      // JOINT LRC-BART MODIFICATION END

      if(update_tau0) {
        const double post_nu=nu0+(double)K0;
        const double post_ssq=(nu0*s0_sq+S0)/post_nu;
        tau0_sq=rsinvchisq_upper(post_nu,post_ssq,tau1_sq,gen);
      }
      if(!w_fixed) {
        w=gen.beta(aw+(double)K0,bw+(double)(Lg-K0));
        w=std::max(1e-12,std::min(1.0-1e-12,w));
      }
    }

    // JOINT LRC-BART ADDITION START
    // PR#2 update_tau/tau_moments: proper normal coefficient conditional,
    // after g hyperparameters and before source variances and AFT augmentation.
    if(joint_model) {
      double precision=1.0/tau_prior_var, weighted_sum=0.0;
      for(size_t i=0;i<n;i++) if(A[i]!=0.0) {
        const double base=fforest.f(i)+((source[i]==1 && Hg>0) ? gforest.f(i) : 0.0);
        precision+=A[i]*A[i]/sigma1_sq;
        weighted_sum+=A[i]*(latent[i]-base)/sigma1_sq;
      }
      const double variance=1.0/precision;
      tau=variance*weighted_sum+std::sqrt(variance)*gen.normal();
    }
    // JOINT LRC-BART ADDITION END

    double rss1=0.0,rss2=0.0;
    for(size_t i=0;i<n;i++) {
      // JOINT LRC-BART MODIFICATION START
      double mu=fforest.f(i)+((source[i]==1 && Hg>0) ? gforest.f(i) : 0.0);
      if(joint_model) mu+=A[i]*tau;
      // JOINT LRC-BART MODIFICATION END
      const double resid=latent[i]-mu;
      if(source[i]==1) rss1+=resid*resid;
      else rss2+=resid*resid;
    }
    if(n2>0) sigma2_sq=(nu_sigma*lambda_sigma2+rss2)/gen.chi_square(nu_sigma+n2);
    if(n1>0) sigma1_sq=(nu_sigma*lambda_sigma1+rss1)/gen.chi_square(nu_sigma+n1);
    else sigma1_sq=sigma2_sq; // approved mBART.real single-arm behavior

    if(augmentation) {
      for(size_t i=0;i<n;i++) if(event[i]==0) {
        // JOINT LRC-BART MODIFICATION START
        // PR#2 impute_censored: treatment also shifts the latent log-time mean.
        double mu=fforest.f(i)+((source[i]==1 && Hg>0) ? gforest.f(i) : 0.0);
        if(joint_model) mu+=A[i]*tau;
        // JOINT LRC-BART MODIFICATION END
        const double sd=std::sqrt(source[i]==1 ? sigma1_sq : sigma2_sq);
        latent[i]=rtnorm(mu,lower[i],sd,gen);
      }
    }

    if(it>=n_burn && ((it-n_burn)%thin==0)) {
      // JOINT LRC-BART ADDITION START
      if(joint_model) out_tau[saved]=tau;
      // JOINT LRC-BART ADDITION END
      out_sig1[saved]=sigma1_sq;
      out_sig2[saved]=sigma2_sq;
      out_tau0[saved]=tau0_sq;
      out_w[saved]=w;
      out_K0[saved]=K0;
      out_Lg[saved]=Lg;

      if(keep_train) {
        for(size_t i=0;i<n;i++) {
          f_train(saved,i)=fforest.f(i);
          g_train(saved,i)=Hg>0 ? gforest.f(i) : 0.0;
        }
      }
      if(keep_test && np>0) {
        fforest.predict(p,np,&x_test[0],&pred_f[0]);
        if(Hg>0) gforest.predict(p,np,&x_test[0],&pred_g[0]);
        for(size_t i=0;i<np;i++) {
          f_test(saved,i)=pred_f[i];
          g_test(saved,i)=Hg>0 ? pred_g[i] : 0.0;
        }
      }
      if(keep_trees) {
        // LRC-BART MODIFICATION: this is the original retained-tree block,
        // extended to save the scalar theta and g-leaf z needed by prediction
        // and realized ESS.  No separate serialization abstraction is used.
        Rcpp::List f_one(Hf);
        for(size_t h=0;h<Hf;h++) {
          tree::cnpv nodes;
          fforest.gettree(h).getnodes(nodes);
          Rcpp::NumericMatrix M(nodes.size(),6);
          for(size_t k=0;k<nodes.size();k++) {
            M(k,0)=(double)nodes[k]->nid();
            M(k,1)=(double)nodes[k]->getv();
            M(k,2)=(double)nodes[k]->getc();
            M(k,3)=nodes[k]->gettheta();
            M(k,4)=(double)nodes[k]->getz();
            M(k,5)=nodes[k]->isleaf() ? 1.0 : 0.0;
          }
          Rcpp::colnames(M)=Rcpp::CharacterVector::create(
            "node_id","variable","cut_index","theta","z","is_leaf");
          f_one[h]=M;
        }
        f_draws[saved]=f_one;

        if(Hg>0) {
          Rcpp::List g_one(Hg);
          for(size_t h=0;h<Hg;h++) {
            tree::cnpv nodes;
            gforest.gettree(h).getnodes(nodes);
            Rcpp::NumericMatrix M(nodes.size(),6);
            for(size_t k=0;k<nodes.size();k++) {
              M(k,0)=(double)nodes[k]->nid();
              M(k,1)=(double)nodes[k]->getv();
              M(k,2)=(double)nodes[k]->getc();
              M(k,3)=nodes[k]->gettheta();
              M(k,4)=(double)nodes[k]->getz();
              M(k,5)=nodes[k]->isleaf() ? 1.0 : 0.0;
            }
            Rcpp::colnames(M)=Rcpp::CharacterVector::create(
              "node_id","variable","cut_index","theta","z","is_leaf");
            g_one[h]=M;
          }
          g_draws[saved]=g_one;
        }
      }
      ++saved;
    }

    if(verbose && ((it+1)%100==0))
      Rcpp::Rcout << "iter " << it+1 << " / " << total << endl;
    if((it%50)==0) Rcpp::checkUserInterrupt();
  }

  Rcpp::NumericMatrix accept(2,3);
  for(size_t k=0;k<3;k++) {
    accept(0,k)=fforest.getaccept(k);
    accept(1,k)=Hg>0 ? gforest.getaccept(k) : NA_REAL;
  }
  Rcpp::rownames(accept)=Rcpp::CharacterVector::create("f","g");
  Rcpp::colnames(accept)=Rcpp::CharacterVector::create("grow","prune","change");

  // Inherited cmbart.cpp return pattern: convert xinfo back to an R list next
  // to the returned tree draws.
  Rcpp::List cutpoints_out(xi.size());
  for(size_t v=0;v<xi.size();v++) cutpoints_out[v]=Rcpp::wrap(xi[v]);

  // JOINT LRC-BART MODIFICATION START
  // Keep the complete legacy output contract when treatment is omitted.
  Rcpp::List result=Rcpp::List::create(
    Rcpp::Named("f_draws")=f_draws,
    Rcpp::Named("g_draws")=g_draws,
    Rcpp::Named("cutpoints")=cutpoints_out,
    Rcpp::Named("sigma1_sq")=out_sig1,
    Rcpp::Named("sigma2_sq")=out_sig2,
    Rcpp::Named("tau0_sq")=out_tau0,
    Rcpp::Named("tau1_sq")=Rcpp::NumericVector(n_draw,tau1_sq),
    Rcpp::Named("w")=out_w,
    Rcpp::Named("K0")=out_K0,
    Rcpp::Named("L_g")=out_Lg,
    Rcpp::Named("f_train")=f_train,
    Rcpp::Named("g_train")=g_train,
    Rcpp::Named("f_test")=f_test,
    Rcpp::Named("g_test")=g_test,
    Rcpp::Named("accept")=accept,
    Rcpp::Named("single_arm")=single_arm);
  if(joint_model) {
    result["tau"]=out_tau; // psi in the manuscript; add A*tau to f+S*g in R
    result["tau_prior_var"]=tau_prior_var;
    result["joint_model"]=true;
    result["no_trial_data"]=no_trial_data;
  }
  return result;
  // JOINT LRC-BART MODIFICATION END
}
// LRC-BART MODIFICATION END

// ---------------------------------------------------------------------------
// LRC-BART ADDITION: prediction from kept pointer-tree serializations.  This
// is a production consumer of theta and z; no test-only leaf exports are kept.
// x_new follows the same p x n convention as clrcbart().
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List clrcbart_predict(Rcpp::List forest_draws,
                            Rcpp::List cutpoints,
                            Rcpp::NumericMatrix x_new)
{
  const size_t p=(size_t)x_new.nrow();
  const size_t n=(size_t)x_new.ncol();
  const size_t nd=(size_t)forest_draws.size();
  // Same inline R-to-xinfo conversion used by clrcbart(), matching the
  // conversion placement in the original cmbart.cpp entry source.
  if((size_t)cutpoints.size()!=p)
    Rcpp::stop("cutpoints must contain one numeric vector per predictor");
  xinfo xi(p);
  for(size_t v=0;v<p;v++) {
    Rcpp::NumericVector cv=cutpoints[v];
    xi[v].assign(cv.begin(),cv.end());
    if(xi[v].empty()) Rcpp::stop("each predictor needs at least one cutpoint");
  }
  Rcpp::NumericMatrix value(nd,n);
  Rcpp::IntegerMatrix nspike(nd,n);

  for(size_t d=0;d<nd;d++) {
    Rcpp::List forest=forest_draws[d];
    for(int h=0;h<forest.size();h++) {
      Rcpp::NumericMatrix M=forest[h];
      std::map<size_t,size_t> row;
      for(int k=0;k<M.nrow();k++) row[(size_t)M(k,0)]=(size_t)k;
      for(size_t i=0;i<n;i++) {
        // LRC-BART ADDITION: traverse the saved mBART node IDs directly here;
        // this is the minimum new behavior required for arbitrary-profile
        // prediction and reached-leaf ESS diagnostics.
        size_t node_id=1;
        std::map<size_t,size_t>::const_iterator it=row.find(node_id);
        if(it==row.end()) Rcpp::stop("invalid serialized tree: root is missing");
        size_t k=it->second;
        while(M(k,5)<0.5) {
          const size_t v=(size_t)M(k,1);
          const size_t c=(size_t)M(k,2);
          if(v>=xi.size() || c>=xi[v].size())
            Rcpp::stop("invalid serialized tree cutpoint index");
          node_id=(x_new(v,i)<xi[v][c]) ? 2*node_id : 2*node_id+1;
          it=row.find(node_id);
          if(it==row.end()) Rcpp::stop("invalid serialized tree: child node is missing");
          k=it->second;
        }
        value(d,i)+=M(k,3);
        nspike(d,i)+=(int)M(k,4);
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("value")=value,
    Rcpp::Named("nspike")=nspike);
}
