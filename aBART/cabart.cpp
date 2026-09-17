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

#include "include.a/rn.h"
#include "include.a/tree.h"
#include "include.a/treefuns.h"
#include "include.a/info.h"
#include "include.a/bartfuns.h"
#include "include.a/bart.h"
#include "include.a/heterbartfuns.h"
#include "include.a/heterbd.h"
#include "include.a/heterbart.h"
#include "include.a/rtnorm.h"

// [[Rcpp::export]]
Rcpp::List cabart(
   SEXP itype,          //1:wbart, 2:pbart, 3:lbart
   SEXP iin,            //number of observations in training data
   SEXP ip,            //dimension of x
   SEXP inp,           //number of observations in test data
   SEXP iix,            //x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP iiy,            //y, train,  nx1
   SEXP idelta,        //censoring indicator
   SEXP iixp,           //x, test, pxnp (transposed so rows are contiguous in memory)
   SEXP im,            //number of trees
   SEXP inc,           //number of cut points
   SEXP ind,           //number of kept draws (except for thinnning ..)
   SEXP iburn,         //number of burn-in draws skipped
   SEXP ithin,         //thinning
   SEXP ipower,
   SEXP ibase,
   SEXP iOffset,
   SEXP itau,
   SEXP inu,
   SEXP ilambda,
   SEXP isigest,
   SEXP iiw,
   SEXP idart,         //dart prior: true(1)=yes, false(0)=no
   SEXP itheta,
   SEXP iomega,
   SEXP igrp,
   SEXP ia,            //param a for sparsity prior
   SEXP ib,            //param b for sparsity prior
   SEXP irho,          //param rho for sparsity prior (default to p)
   SEXP iaug,          //categorical strategy: true(1)=data augment false(0)=degenerate trees
   SEXP inprintevery,
   SEXP iXinfo,
   SEXP iseed          //seed for random number generation
)
{
   //process args
   int type = Rcpp::as<int>(itype);
   size_t n = Rcpp::as<int>(iin);
   size_t p = Rcpp::as<int>(ip);
   size_t np = Rcpp::as<int>(inp);
   Rcpp::NumericVector  xv(iix);
   double *ix = &xv[0];
   Rcpp::NumericVector  yv(iiy); 
   double *iy = &yv[0];
   Rcpp::IntegerVector  deltav(idelta); 
   int *delta = &deltav[0];
   
   Rcpp::NumericMatrix xpv(iixp);
   double *ixp = nullptr;
   if(np>0) ixp = &xpv[0];
   size_t m = Rcpp::as<int>(im);
   Rcpp::IntegerVector _nc(inc);
   int *numcut = &_nc[0];
   size_t nd = Rcpp::as<int>(ind);
   int burn = Rcpp::as<int>(iburn);
   size_t thin = Rcpp::as<int>(ithin);
   double mybeta = Rcpp::as<double>(ipower);
   double alpha = Rcpp::as<double>(ibase);
   double Offset = Rcpp::as<double>(iOffset);
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
   int seed = Rcpp::as<int>(iseed);

   size_t nkeeptrain = nd/thin;     //Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = nd/thin;      //Rcpp::as<int>(_inkeeptest);
   size_t nkeeptreedraws = nd/thin; //Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(inprintevery);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
   Rcpp::NumericMatrix Xinfo(iXinfo);
   Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericVector acc(nd+burn);
   Rcpp::NumericVector nodes(nd+burn);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);

   //random number generation
   arn gen;
   gen.set_seed(seed);

   heterbart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       for(int j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
       //for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;

   size_t skiptr=thin, skipte=thin, skiptreedraws=thin;
   //--------------------------------------------------
   //create temporaries
   double df=n+nu;
   double *z = new double[n]; 
   double *svec = new double[n]; 
   double *sign = NULL;
   if(type!=1) sign = new double[n]; 

   for(size_t i=0; i<n; i++) {
     if(type==1) {
       svec[i] = iw[i]*sigma; 
       z[i]=iy[i]; 
     }
     else {
       svec[i] = 1.;
       if(iy[i]==0) sign[i] = -1.;
       else sign[i] = 1.;
       z[i] = sign[i];
     }
   }
   //--------------------------------------------------
   //set up BART model
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,z,numcut);

   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);
   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; 
   if(np) { fhattest = new double[np]; }

   //--------------------------------------------------
   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
   bool keeptest,/*keeptestme*/keeptreedraw;

   time_t tp;
   int time1 = time(&tp), total=nd+burn;
   xinfo& xi = bm.getxinfo();
   std::vector<double> allLsum(m);
   
   for(int i=0;i<total;i++) {
     // cout << "--- Iteration " << i << " ---"<< endl;
     fill(allLsum.begin(), allLsum.end(), 0.0);
     
      bm.draw(svec,allLsum,gen);
      double L = accumulate(allLsum.begin(), 
                            allLsum.end(), 
                            0.0);
      
      if(type==1) {
        //draw sigma
        double rss=0.;
        for(size_t k=0;k<n;k++) rss += pow((z[k]-bm.f(k))/(iw[k]), 2.); 
        sigma = sqrt((nu*lambda + rss)/gen.chi_square(df));
        sdraw[i]=sigma;
      }

      for(size_t k=0; k<n; k++) {
        if(type==1) {
          svec[k]=iw[k]*sigma;
          if(delta[k]==0) z[k]= rtnorm(bm.f(k), iy[k], svec[k], gen);
        }
      }
      
      acc[i] = bm.getaccept();
      nodes[i] = L;
      if(i>=burn) {
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=Offset+bm.f(k);
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         if(keeptest) {
	   bm.predict(p,np,ixp,fhattest);
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=Offset+fhattest[k];
            tecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);

	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
	    size_t k=(i-burn)/skiptreedraws;
	    for(size_t j=0;j<p;j++){
	      varcnt(k,j)=ivarcnt[j];
	      varprb(k,j)=ivarprb[j];
	    }
	    }
         }
      }
   }
   int time2 = time(&tp);

   if(fhattest) delete[] fhattest;
   delete[] z;
   delete[] svec;
   if(type!=1) delete[] sign;

   //return list
   Rcpp::List ret;
   if(type==1) ret["sigma"]=sdraw;
   ret["yhat.train"]=trdraw;
   ret["yhat.test"]=tedraw;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;
   ret["accept"]=acc;
   ret["L"]=nodes;
   
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
