// cess.cpp
// Centralized compiled ESS source for LRC-BART.
//
// BASELINE AND PROVENANCE
// This source takes the sourceCpp entry style, pointer-tree traversal, and
// per-leaf profile counting from the original project files
//   mapbart-sim-gaussian/ess_local/cwbart_ess.cpp
//   mapbart-sim-survival/ess_local/cabart_ess.cpp
// and changes only the compiled responsibility required by the revision:
// prior-tree c_g, ESS_tau0 evaluation, and exact posterior partition variance.
// Stage-1 fitting, the Pr-rule, study configuration, and plotting remain R-side
// responsibilities of each subproject's inherited ess_local workflow.

#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

typedef std::vector<double> v1d;

#include "include.lrc/rn.h"
#include "include.lrc/tree.h"
#include "include.lrc/treefuns.h"

// ---------------------------------------------------------------------------
// LRC-BART ADDITION START
// Helpers use the same xinfo, pointer tree, ancestor bounds, and node-id
// conventions as clrcbart.cpp.
// ---------------------------------------------------------------------------
static xinfo cess_xinfo_from_list(const Rcpp::List& cutpoints, size_t p)
{
  if((size_t)cutpoints.size()!=p)
    Rcpp::stop("cutpoints must contain one numeric vector per predictor");
  xinfo xi(p);
  for(size_t v=0;v<p;v++) {
    Rcpp::NumericVector cv=cutpoints[v];
    xi[v].assign(cv.begin(),cv.end());
    if(xi[v].empty()) Rcpp::stop("each predictor needs at least one cutpoint");
  }
  return xi;
}

static void cess_grow_prior(tree& t, tree::tree_p node, size_t depth,
                            double alpha, double beta,
                            xinfo& xi, size_t maxdepth, rn& gen)
{
  if(depth>=maxdepth) return;
  if(gen.uniform()>=alpha/std::pow(1.0+(double)depth,beta)) return;
  std::vector<size_t> goodvars;
  getgoodvars(node,xi,goodvars);
  if(goodvars.empty()) return;
  const size_t v=goodvars[(size_t)std::floor(gen.uniform()*goodvars.size())];
  int L=0;
  int U=(int)xi[v].size()-1;
  node->rg(v,&L,&U);
  const size_t c=(size_t)(L+(int)std::floor(gen.uniform()*(U-L+1)));
  t.birthp(node,v,c,0.0,0.0,1,1);
  cess_grow_prior(t,node->getl(),depth+1,alpha,beta,xi,maxdepth,gen);
  cess_grow_prior(t,node->getr(),depth+1,alpha,beta,xi,maxdepth,gen);
}

static size_t cess_serialized_row(const std::map<size_t,size_t>& row,
                                  size_t node_id)
{
  std::map<size_t,size_t>::const_iterator it=row.find(node_id);
  if(it==row.end()) Rcpp::stop("invalid serialized tree: child node is missing");
  return it->second;
}

static size_t cess_serialized_leaf(const Rcpp::NumericMatrix& M,
                                   const std::map<size_t,size_t>& row,
                                   const xinfo& xi,
                                   const Rcpp::NumericMatrix& x,
                                   size_t i)
{
  size_t node_id=1;
  size_t k=cess_serialized_row(row,node_id);
  while(M(k,5)<0.5) {
    const size_t v=(size_t)M(k,1);
    const size_t c=(size_t)M(k,2);
    if(v>=xi.size() || c>=xi[v].size())
      Rcpp::stop("invalid serialized tree cutpoint index");
    node_id=(x(v,i)<xi[v][c]) ? 2*node_id : 2*node_id+1;
    k=cess_serialized_row(row,node_id);
  }
  return k;
}

static double cess_value(const Rcpp::NumericVector& x, size_t i, size_t n)
{
  if(x.size()==1) return x[0];
  if((size_t)x.size()!=n) Rcpp::stop("ESS vector length must be one or the number of draws");
  return x[i];
}
// LRC-BART ADDITION END

// ---------------------------------------------------------------------------
// LRC-BART ADDITION: c_g = E_T sum_l (n_l/N)^2 under one prior g tree.
// x follows mBART's p x n convention.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
double cess_cg(Rcpp::NumericMatrix x, Rcpp::List cutpoints,
               double alpha_g=0.5, double beta_g=3.0,
               int n_sim=2000, int max_depth=12, int seed=1)
{
  const size_t p=(size_t)x.nrow();
  const size_t n=(size_t)x.ncol();
  if(p==0 || n==0 || n_sim<1 || max_depth<1) Rcpp::stop("invalid c_g inputs");
  xinfo xi=cess_xinfo_from_list(cutpoints,p);
  arn gen;
  gen.set_seed((unsigned int)seed);
  double total=0.0;

  for(int s=0;s<n_sim;s++) {
    tree t;
    cess_grow_prior(t,&t,0,alpha_g,beta_g,xi,(size_t)max_depth,gen);
    tree::npv leaves;
    t.getbots(leaves);
    std::map<tree::tree_p,size_t> index;
    for(size_t l=0;l<leaves.size();l++) index[leaves[l]]=l;
    std::vector<double> count(leaves.size(),0.0);
    for(size_t i=0;i<n;i++) count[index[t.bn(&x[0]+i*p,xi)]]+=1.0;
    double one=0.0;
    for(size_t l=0;l<count.size();l++) {
      const double share=count[l]/(double)n;
      one+=share*share;
    }
    total+=one;
    if((s%100)==0) Rcpp::checkUserInterrupt();
  }
  return total/(double)n_sim;
}

// ---------------------------------------------------------------------------
// LRC-BART ADDITION: evaluate the Stage-1 ESS_tau0 formula for each V_mu^f block.
// The expectation over tau0^2 is the mean over tau0_sq.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector cess_ess_tau0(Rcpp::NumericVector V_mu_f,
                                  Rcpp::NumericVector sigma1_sq,
                                  double c_g, int H_g,
                                  Rcpp::NumericVector tau0_sq)
{
  if(V_mu_f.size()==0 || tau0_sq.size()==0 || H_g<0 || c_g<0.0)
    Rcpp::stop("invalid ESS_tau0 inputs");
  const size_t nb=(size_t)V_mu_f.size();
  if(sigma1_sq.size()!=1 && (size_t)sigma1_sq.size()!=nb)
    Rcpp::stop("sigma1_sq must have length one or length(V_mu_f)");
  Rcpp::NumericVector out(nb);
  for(size_t b=0;b<nb;b++) {
    const double sig=sigma1_sq.size()==1 ? sigma1_sq[0] : sigma1_sq[b];
    double mean_inv=0.0;
    for(int j=0;j<tau0_sq.size();j++)
      mean_inv+=1.0/(V_mu_f[b]+c_g*(double)H_g*tau0_sq[j]);
    out[b]=sig*mean_inv/(double)tau0_sq.size();
  }
  return out;
}

// ---------------------------------------------------------------------------
// LRC-BART ADDITION: exact realized discrepancy variance
//   sum_h sum_l (n_l/N)^2 tau^2_{z_hl}
// for every kept posterior g-forest draw.  This is the compiled quantity used
// by scalar, regional, and profile-subset ESS reporting.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector cess_forest_gvar(Rcpp::List g_draws,
                                     Rcpp::List cutpoints,
                                     Rcpp::NumericMatrix x_profiles,
                                     Rcpp::NumericVector tau0_sq,
                                     Rcpp::NumericVector tau1_sq)
{
  const size_t nd=(size_t)g_draws.size();
  const size_t p=(size_t)x_profiles.nrow();
  const size_t n=(size_t)x_profiles.ncol();
  if(nd==0 || p==0 || n==0) Rcpp::stop("g draws and profiles must be nonempty");
  xinfo xi=cess_xinfo_from_list(cutpoints,p);
  Rcpp::NumericVector out(nd);

  for(size_t d=0;d<nd;d++) {
    Rcpp::List forest=g_draws[d];
    double Vg=0.0;
    for(int h=0;h<forest.size();h++) {
      Rcpp::NumericMatrix M=forest[h];
      std::map<size_t,size_t> row;
      for(int k=0;k<M.nrow();k++) row[(size_t)M(k,0)]=(size_t)k;
      std::vector<double> count(M.nrow(),0.0);
      for(size_t i=0;i<n;i++)
        count[cess_serialized_leaf(M,row,xi,x_profiles,i)]+=1.0;
      for(int k=0;k<M.nrow();k++) {
        if(M(k,5)<0.5 || count[k]==0.0) continue;
        const double share=count[k]/(double)n;
        const double leaf_variance=M(k,4)>0.5 ?
          cess_value(tau0_sq,d,nd) : cess_value(tau1_sq,d,nd);
        Vg+=share*share*leaf_variance;
      }
    }
    out[d]=Vg;
  }
  return out;
}
