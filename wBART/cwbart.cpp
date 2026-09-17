#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>

// Create aliases 
#define printf Rprintf
#define cout Rcpp::Rcout
#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

typedef std::vector<double> v1d; //1D vector
typedef std::vector<v1d> v2d; //2D vector (2x2 matrix)
typedef std::vector<v2d> v3d; //2x2x2 3D vector

using std::endl;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "include.w/rn.h"
#include "include.w/tree.h"
#include "include.w/treefuns.h"
#include "include.w/info.h"
#include "include.w/bartfuns.h"
#include "include.w/bd.h"
#include "include.w/bart.h"
#include "include.w/heterbartfuns.h"
#include "include.w/heterbd.h"
#include "include.w/heterbart.h"

// [[Rcpp::export]]
Rcpp::List cwbart(
    SEXP iin,            //number of observations in training data
    SEXP ip,		//dimension of x
    SEXP inp,		//number of observations in test data
    SEXP iix,		//x, train,  pxn (transposed so rows are contiguous in memory)
    SEXP iiy,		//y, train,  nx1
    SEXP iixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
    SEXP im,		//number of trees
    SEXP inc,		//number of cut points
    SEXP ind,		//number of kept draws (except for thinnning ..)
    SEXP iburn,		//number of burn-in draws skipped
    SEXP ipower,
    SEXP ibase,
    SEXP itau,
    SEXP inu,
    SEXP ilambda,
    SEXP isigest,
    SEXP iiw,
    SEXP idart,
    SEXP itheta,
    SEXP iomega,
    SEXP igrp,
    SEXP ia,
    SEXP ib,
    SEXP irho,
    SEXP iaug,
    SEXP inkeeptrain,
    SEXP inkeeptest,
    SEXP inkeeptestme,
    SEXP inkeeptreedraws,
    SEXP inprintevery,
    SEXP iXinfo
)
{
  size_t n = Rcpp::as<int>(iin);
  size_t p = Rcpp::as<int>(ip);
  size_t np = Rcpp::as<int>(inp);
  Rcpp::NumericVector  xv(iix);
  double *ix = &xv[0];
  Rcpp::NumericVector  yv(iiy); 
  double *iy = &yv[0];
  Rcpp::NumericMatrix xpv(iixp);
  double *ixp = nullptr;
  if(np>0) ixp = &xpv[0];
  size_t m = Rcpp::as<int>(im);
  Rcpp::IntegerVector _nc(inc);
  int *numcut = &_nc[0];
  size_t nd = Rcpp::as<int>(ind);
  size_t burn = Rcpp::as<int>(iburn);
  double mybeta = Rcpp::as<double>(ipower);
  double alpha = Rcpp::as<double>(ibase);
  double tau = Rcpp::as<double>(itau);
  double nu = Rcpp::as<double>(inu);
  double lambda = Rcpp::as<double>(ilambda);
  double sigma=Rcpp::as<double>(isigest);
  Rcpp::NumericVector  wv(iiw); 
  double *iw = &wv[0];
  bool dart;
  if(Rcpp::as<int>(idart)==1) dart=true;
  else dart=false;
  double a = Rcpp::as<double>(ia);
  double b = Rcpp::as<double>(ib);
  double rho = Rcpp::as<double>(irho);
  bool aug;
  if(Rcpp::as<int>(iaug)==1) aug=true;
  else aug=false;
  double theta = Rcpp::as<double>(itheta);
  double omega = Rcpp::as<double>(iomega);
  Rcpp::IntegerVector _grp(igrp);
  size_t nkeeptrain = Rcpp::as<int>(inkeeptrain);
  size_t nkeeptest = Rcpp::as<int>(inkeeptest);
  size_t nkeeptestme = Rcpp::as<int>(inkeeptestme);
  size_t nkeeptreedraws = Rcpp::as<int>(inkeeptreedraws);
  size_t printevery = Rcpp::as<int>(inprintevery);
  Rcpp::NumericMatrix Xinfo(iXinfo);
  
  Rcpp::NumericVector trmean(n); //train
  Rcpp::NumericVector temean(np);
  Rcpp::NumericVector sdraw(nd+burn);
  Rcpp::NumericMatrix trdraw(nkeeptrain,n);
  Rcpp::NumericMatrix tedraw(nkeeptest,np);
  Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
  Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
  
  Rcpp::NumericVector acc(nd+burn);
  Rcpp::NumericMatrix node(nd+burn,m);
  Rcpp::NumericVector nodesum(nd+burn);
  
  arn gen;
  
  heterbart bm(m);
  
  if(Xinfo.size()>0) {
    xinfo _xi;
    _xi.resize(p);
    for(size_t i=0;i<p;i++) {
      _xi[i].resize(numcut[i]);
      for(int j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
    }
    bm.setxinfo(_xi);
  }

    for(size_t i=0;i<n;i++) trmean[i]=0.0;
    for(size_t i=0;i<np;i++) temean[i]=0.0;
    
    //-----------------------------------------------------------
    
    size_t skiptr,skipte,skipteme,skiptreedraws;
    if(nkeeptrain) {skiptr=nd/nkeeptrain;}
    else skiptr = nd+1;
    if(nkeeptest) {skipte=nd/nkeeptest;}
    else skipte=nd+1;
    if(nkeeptestme) {skipteme=nd/nkeeptestme;}
    else skipteme=nd+1;
    if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
    else skiptreedraws=nd+1;
    
    //--------------------------------------------------
    bm.setprior(alpha,mybeta,tau);
    bm.setdata(p,n,ix,iy,numcut);
    bm.setdart(a,b,rho,aug,dart,theta,omega);
    
    //--------------------------------------------------
    double *svec = new double[n];
    for(size_t i=0;i<n;i++) svec[i]=iw[i]*sigma;
    
    //--------------------------------------------------
    std::stringstream treess;  //string stream to write trees to  
    treess.precision(10);
    treess << nkeeptreedraws << " " << m << " " << p << endl;
    // dart iterations
    std::vector<double> ivarprb (p,0.);
    std::vector<size_t> ivarcnt (p,0);
    
    //--------------------------------------------------
    double* fhattest=0; //posterior mean for prediction
    if(np) { fhattest = new double[np]; }
    double restemp=0.0,rss=0.0;
    
    //--------------------------------------------------
    size_t trcnt=0; //count kept train draws
    size_t tecnt=0; //count kept test draws
    size_t temecnt=0; //count test draws into posterior mean
    size_t treedrawscnt=0; //count kept bart draws
    bool keeptest,keeptestme,keeptreedraw;
    
    xinfo& xi = bm.getxinfo();
    
    std::vector<double> allLsum(m);
    
    for(size_t i=0;i<(nd+burn);i++) {
      
      if(i==(burn/2)&&dart) bm.startdart();
      
      // cout << "--- Iteration " << i << " ---" << endl;
      
      fill(allLsum.begin(), allLsum.end(), 0.0);
      bm.draw(svec,
              allLsum,
              gen);
      
      acc[i] = bm.getaccept();
      double L = accumulate(allLsum.begin(),
                            allLsum.end(),
                            0.0);
      for(size_t j=0;j<m;j++) {
        node(i,j) = allLsum[j];
      }
      nodesum[i] = L;
      
      rss=0.0;
      for(size_t k=0;k<n;k++) {restemp=(iy[k]-bm.f(k))/(iw[k]); rss += restemp*restemp;}
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(n+nu));
      for(size_t k=0;k<n;k++) svec[k]=iw[k]*sigma;
      sdraw[i]=sigma;
      if(i>=burn) {
        for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
        if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
          for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
          trcnt+=1;
        }
        keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
        keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
        if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
        if(keeptest) {
          for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
          tecnt+=1;
        }
        if(keeptestme) {
          for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
          temecnt+=1;
        }
        keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
        if(keeptreedraw) {
          for(size_t j=0;j<m;j++) {
            treess << bm.gettree(j);
          }
          ivarcnt=bm.getnv();
          ivarprb=bm.getpv();
          size_t k=(i-burn)/skiptreedraws;
          for(size_t j=0;j<p;j++){
            varcnt(k,j)=ivarcnt[j];
            varprb(k,j)=ivarprb[j];
          }
          treedrawscnt +=1;
        }
      }
    }
    for(size_t k=0;k<n;k++) trmean[k]/=nd;
    for(size_t k=0;k<np;k++) temean[k]/=temecnt;

    if(fhattest) delete[] fhattest;
    if(svec) delete [] svec;
    
    //--------------------------------------------------
    Rcpp::List ret;
    ret["sigma"]=sdraw;
    ret["yhat.train.mean"]=trmean;
    ret["yhat.train"]=trdraw;
    ret["yhat.test.mean"]=temean;
    ret["yhat.test"]=tedraw;
    ret["varcount"]=varcnt;
    ret["varprob"]=varprb;
    
    ret["accept"]=acc;
    ret["node"]=node;
    ret["L"]=nodesum;

    Rcpp::List xiret(xi.size());
    for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
    }
    
    Rcpp::List treesL;
    treesL["cutpoints"] = xiret;
    treesL["trees"]=Rcpp::CharacterVector(treess.str());
    ret["treedraws"] = treesL;
    
    return ret;
  }